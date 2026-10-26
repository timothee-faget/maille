use crate::assembly::{
    element_conduction, element_dofs, element_heat_fluxes, element_heat_source, element_stiffness,
    element_stresses, face_convection, FaceConvection,
};
use crate::element::{bar_equivalent_load, beam_equivalent_load};
use crate::error::MailleError;
use crate::material::positive;
use crate::mesh::{Dof, Element, Mesh};
use crate::results::Field;
use crate::solve::{reverse_cuthill_mckee, Cholesky, SkylineMatrix};

/// Condition aux limites de Dirichlet : valeur imposée à un DOF d'un Node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bc {
    pub node: usize,
    pub dof: Dof,
    pub value: f64,
}

/// Force (ou moment pour `Theta`) nodale, Neumann. En thermique (`Dof::T`) : puissance apportée au
/// Node (négative si elle sort).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Load {
    pub node: usize,
    pub dof: Dof,
    pub value: f64,
}

/// Charge répartie uniforme sur un Element, par unité de longueur (Bar : axiale ; Beam : dans le sens de V).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DistributedLoad {
    pub element: usize,
    pub value: f64,
}

/// Source de chaleur uniforme dans un Element (thermique) : puissance par unité de volume
/// (ex. W/mm³ avec des W et des mm).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct HeatSource {
    pub element: usize,
    pub value: f64,
}

/// Convection h (T − T∞) sur une face d'un Element (thermique). Numérotation des faces :
/// voir `Element::face_nodes`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Convection {
    pub element: usize,
    pub face: usize,
    /// Coefficient d'échange h, strictement positif : puissance par unité de surface et par degré
    /// (ex. W/(mm²·K)). Nom Abaqus : *film*.
    pub film: f64,
    /// Température du fluide T∞, quelconque mais finie. Nom Abaqus : *sink*.
    pub sink: f64,
}

/// Physique résolue sur le Mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Analysis {
    /// Statique linéaire : déplacements, contraintes.
    Static,
    /// Thermique stationnaire : −∇·(k ∇T) = Q, un DOF `T` par Node.
    Thermal,
}

#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    mesh: Mesh,
    analysis: Analysis,
    bcs: Vec<Bc>,
    loads: Vec<Load>,
    distributed: Vec<DistributedLoad>,
    heat_sources: Vec<HeatSource>,
    convection: Vec<Convection>,
}

impl Model {
    /// Analyse statique linéaire, voir `with_analysis`.
    pub fn new(mesh: Mesh, bcs: Vec<Bc>, loads: Vec<Load>) -> Result<Self, MailleError> {
        Self::with_analysis(mesh, Analysis::Static, bcs, loads)
    }

    /// Le Material de chaque Element doit avoir les propriétés de l'analyse : en statique E, et ν
    /// pour les triangles et tétraèdres ; en thermique k (triangles et tétraèdres seulement).
    pub fn with_analysis(
        mesh: Mesh,
        analysis: Analysis,
        bcs: Vec<Bc>,
        loads: Vec<Load>,
    ) -> Result<Self, MailleError> {
        for (i, e) in mesh.elements.iter().enumerate() {
            let m = e.material();
            let one_d = matches!(e, Element::Bar(_) | Element::Beam(_));
            match analysis {
                Analysis::Static if m.young.is_none() => {
                    return Err(MailleError::MissingYoung { element: i })
                }
                Analysis::Static if !one_d && m.poisson.is_none() => {
                    return Err(MailleError::MissingPoisson { element: i })
                }
                Analysis::Thermal if one_d => {
                    return Err(MailleError::UnsupportedThermalElement { element: i })
                }
                Analysis::Thermal if m.conductivity.is_none() => {
                    return Err(MailleError::MissingConductivity { element: i })
                }
                _ => {}
            }
        }
        let mut seen = vec![false; mesh.nodes.len() * dofs_per_node(&mesh, analysis)];
        for bc in &bcs {
            let i = dof_index(&mesh, analysis, bc.node, bc.dof)?;
            if std::mem::replace(&mut seen[i], true) {
                return Err(MailleError::DuplicateBc {
                    node: bc.node,
                    dof: bc.dof,
                });
            }
        }
        for l in &loads {
            dof_index(&mesh, analysis, l.node, l.dof)?;
        }
        Ok(Self {
            mesh,
            analysis,
            bcs,
            loads,
            distributed: Vec::new(),
            heat_sources: Vec::new(),
            convection: Vec::new(),
        })
    }

    pub fn with_distributed_loads(
        mut self,
        distributed: Vec<DistributedLoad>,
    ) -> Result<Self, MailleError> {
        for d in &distributed {
            self.element(d.element)?;
        }
        self.distributed = distributed;
        Ok(self)
    }

    /// Sources de chaleur volumiques, en analyse thermique seulement.
    pub fn with_heat_sources(mut self, sources: Vec<HeatSource>) -> Result<Self, MailleError> {
        self.require_thermal()?;
        for s in &sources {
            self.element(s.element)?;
        }
        self.heat_sources = sources;
        Ok(self)
    }

    /// Convection sur des faces, en analyse thermique seulement.
    pub fn with_convection(mut self, convection: Vec<Convection>) -> Result<Self, MailleError> {
        self.require_thermal()?;
        for c in &convection {
            if self.element(c.element)?.face_positions(c.face).is_none() {
                return Err(MailleError::UnknownFace {
                    element: c.element,
                    face: c.face,
                });
            }
            positive("h (coefficient d'échange)", c.film)?;
            if !c.sink.is_finite() {
                return Err(MailleError::NonFinite {
                    name: "T∞ (température du fluide)",
                    value: c.sink,
                });
            }
        }
        self.convection = convection;
        Ok(self)
    }

    /// L'Element `i`, ou `UnknownElement`.
    fn element(&self, i: usize) -> Result<&Element, MailleError> {
        self.mesh
            .elements
            .get(i)
            .ok_or(MailleError::UnknownElement { element: i })
    }

    fn require_thermal(&self) -> Result<(), MailleError> {
        if self.analysis == Analysis::Thermal {
            Ok(())
        } else {
            Err(MailleError::ThermalOnly)
        }
    }

    /// Résout K_ff u_f = F_f − K_fc u_c (K : rigidité en statique, conduction en thermique), puis
    /// réactions R = K u − F sur les DOF contraints. K n'est jamais formée en entier : K_ff est
    /// assemblée en ligne de ciel, les réactions sont recalculées Element par Element.
    pub fn solve(&self) -> Result<Field, MailleError> {
        let mesh = &self.mesh;
        let n = mesh.nodes.len() * dofs_per_node(mesh, self.analysis);
        let dofs: Vec<Vec<usize>> = mesh
            .elements
            .iter()
            .map(|e| match self.analysis {
                Analysis::Static => element_dofs(e),
                Analysis::Thermal => e.nodes().to_vec(),
            })
            .collect();
        // convection : termes de face, regroupés par Element, ajoutés à sa matrice et à F
        let mut faces: Vec<Vec<FaceConvection>> = vec![Vec::new(); mesh.elements.len()];
        for c in &self.convection {
            if let Some(terms) = face_convection(mesh, c.element, c.face, c.film, c.sink) {
                faces[c.element].push(terms);
            }
        }
        let matrix = |i: usize| match self.analysis {
            Analysis::Static => element_stiffness(mesh, i),
            Analysis::Thermal => {
                let mut k = element_conduction(mesh, i);
                for face in &faces[i] {
                    face.add_matrix(&mut k);
                }
                k
            }
        };

        let mut f = vec![0.0; n];
        for l in &self.loads {
            f[dof_index(mesh, self.analysis, l.node, l.dof)?] += l.value;
        }
        for d in &self.distributed {
            let e = &mesh.elements[d.element];
            let equivalent = match e {
                Element::Bar(_) => bar_equivalent_load(d.value, mesh.length(d.element)).to_vec(),
                Element::Beam(_) => beam_equivalent_load(d.value, mesh.length(d.element)).to_vec(),
                Element::Tri3(_) | Element::Tri6(_) | Element::Tet4(_) | Element::Tet10(_) => {
                    return Err(MailleError::UnsupportedDistributedLoad { element: d.element })
                }
            };
            for (&g, v) in dofs[d.element].iter().zip(equivalent) {
                f[g] += v;
            }
        }
        for s in &self.heat_sources {
            let equivalent = element_heat_source(mesh, s.element, s.value);
            for (&g, v) in dofs[s.element].iter().zip(equivalent) {
                f[g] += v;
            }
        }
        for (terms, element_dofs) in faces.iter().zip(&dofs) {
            for face in terms {
                face.add_load(element_dofs, &mut f);
            }
        }
        let mut u = vec![0.0; n];
        let mut fixed = vec![false; n];
        for bc in &self.bcs {
            let i = dof_index(mesh, self.analysis, bc.node, bc.dof)?;
            u[i] = bc.value;
            fixed[i] = true;
        }

        FreeSystem::new(&dofs, &matrix, &fixed, &u)?.solve(&f, &mut u);

        let mut reactions = vec![0.0; n];
        for (i, dofs) in dofs.iter().enumerate() {
            if !dofs.iter().any(|&g| fixed[g]) {
                continue;
            }
            for (row, &g) in matrix(i).iter().zip(dofs).filter(|(_, &g)| fixed[g]) {
                reactions[g] += row.iter().zip(dofs).map(|(k, &h)| k * u[h]).sum::<f64>();
            }
        }
        for i in (0..n).filter(|&i| fixed[i]) {
            reactions[i] -= f[i];
        }
        let elements = 0..mesh.elements.len();
        Ok(match self.analysis {
            Analysis::Static => Field {
                stresses: elements.map(|i| element_stresses(mesh, i, &u)).collect(),
                displacements: u,
                reactions,
                temperatures: Vec::new(),
                heat_fluxes: Vec::new(),
            },
            Analysis::Thermal => Field {
                heat_fluxes: elements.map(|i| element_heat_fluxes(mesh, i, &u)).collect(),
                temperatures: u,
                reactions,
                displacements: Vec::new(),
                stresses: Vec::new(),
            },
        })
    }
}

/// 1 DOF (T) par Node en thermique, sinon selon la famille d'Elements.
fn dofs_per_node(mesh: &Mesh, analysis: Analysis) -> usize {
    match analysis {
        Analysis::Static => mesh.dofs_per_node(),
        Analysis::Thermal => 1,
    }
}

/// Indice global du DOF `dof` du Node `node` pour cette analyse.
fn dof_index(mesh: &Mesh, analysis: Analysis, node: usize, dof: Dof) -> Result<usize, MailleError> {
    match analysis {
        Analysis::Static => mesh.dof_index(node, dof),
        Analysis::Thermal if node >= mesh.nodes.len() => Err(MailleError::UnknownNode { node }),
        Analysis::Thermal if dof == Dof::T => Ok(node),
        Analysis::Thermal => Err(MailleError::InvalidDof { node, dof }),
    }
}

/// Système K_ff u_f = F_f − K_fc u_c factorisé une fois, pour résoudre autant de seconds membres F
/// qu'on veut tant que K et les valeurs imposées u_c ne changent pas.
struct FreeSystem {
    /// free[g] : indice du DOF global g parmi les DOF libres
    free: Vec<Option<usize>>,
    /// rank[a] : position du DOF libre a dans la numérotation du solveur
    rank: Vec<usize>,
    factor: Cholesky,
    /// K_fc u_c, dans la numérotation du solveur
    coupling: Vec<f64>,
}

impl FreeSystem {
    /// `dofs[i]` : DOF globaux de l'Element i, `matrix(i)` : sa matrice, `u` : valeurs imposées
    /// aux DOF `fixed`. Les DOF libres sont renumérotés par Cuthill-McKee inverse (graphe : deux DOF
    /// sont voisins s'ils partagent un Element) avant l'assemblage en ligne de ciel.
    fn new(
        dofs: &[Vec<usize>],
        matrix: &impl Fn(usize) -> Vec<Vec<f64>>,
        fixed: &[bool],
        u: &[f64],
    ) -> Result<Self, MailleError> {
        let mut free = vec![None; fixed.len()];
        let mut count = 0;
        for (g, _) in fixed.iter().enumerate().filter(|(_, &c)| !c) {
            free[g] = Some(count);
            count += 1;
        }
        let element_free: Vec<Vec<usize>> = dofs
            .iter()
            .map(|d| d.iter().filter_map(|&g| free[g]).collect())
            .collect();

        let mut adjacency = vec![Vec::new(); count];
        for local in &element_free {
            for &a in local {
                adjacency[a].extend(local.iter().filter(|&&b| b != a));
            }
        }
        for neighbours in &mut adjacency {
            neighbours.sort_unstable();
            neighbours.dedup();
        }
        let mut rank = vec![0; count];
        for (r, a) in reverse_cuthill_mckee(&adjacency).into_iter().enumerate() {
            rank[a] = r;
        }
        let mut first: Vec<usize> = (0..count).collect();
        for (a, neighbours) in adjacency.iter().enumerate() {
            for &b in neighbours {
                first[rank[a]] = first[rank[a]].min(rank[b]);
            }
        }

        let mut kff = SkylineMatrix::zeros(first);
        let mut coupling = vec![0.0; count];
        for (i, dofs) in dofs.iter().enumerate() {
            for (row, &gi) in matrix(i).iter().zip(dofs) {
                let Some(a) = free[gi] else { continue };
                for (k, &gj) in row.iter().zip(dofs) {
                    match free[gj] {
                        // K symétrique : chaque paire (i, j) est ajoutée une fois, en dessous de la diagonale
                        Some(b) if rank[b] <= rank[a] => kff.add(rank[a], rank[b], *k),
                        Some(_) => {}
                        None => coupling[rank[a]] += k * u[gj],
                    }
                }
            }
        }
        Ok(Self {
            free,
            rank,
            factor: kff.factorize()?,
            coupling,
        })
    }

    /// Résout K_ff u_f = F_f − K_fc u_c (`f` : vecteur global) et écrit u_f dans `u`.
    fn solve(&self, f: &[f64], u: &mut [f64]) {
        let mut rhs: Vec<f64> = self.coupling.iter().map(|c| -c).collect();
        for (g, a) in self.free.iter().enumerate() {
            if let Some(a) = a {
                rhs[self.rank[*a]] += f[g];
            }
        }
        let x = self.factor.solve(&rhs);
        for (g, a) in self.free.iter().enumerate() {
            if let Some(a) = a {
                u[g] = x[self.rank[*a]];
            }
        }
    }
}
