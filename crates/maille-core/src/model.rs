use crate::assembly::{element_dofs, element_stiffness, element_stresses};
use crate::element::{bar_equivalent_load, beam_equivalent_load};
use crate::error::MailleError;
use crate::mesh::{Dof, Element, Mesh};
use crate::results::Field;
use crate::solve::{reverse_cuthill_mckee, SkylineMatrix};

/// Condition aux limites de Dirichlet : valeur imposée à un DOF d'un Node.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bc {
    pub node: usize,
    pub dof: Dof,
    pub value: f64,
}

/// Force (ou moment pour `Theta`) nodale, Neumann.
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

#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    mesh: Mesh,
    bcs: Vec<Bc>,
    loads: Vec<Load>,
    distributed: Vec<DistributedLoad>,
}

impl Model {
    pub fn new(mesh: Mesh, bcs: Vec<Bc>, loads: Vec<Load>) -> Result<Self, MailleError> {
        let mut seen = vec![false; mesh.nodes.len() * mesh.dofs_per_node()];
        for bc in &bcs {
            let i = mesh.dof_index(bc.node, bc.dof)?;
            if std::mem::replace(&mut seen[i], true) {
                return Err(MailleError::DuplicateBc {
                    node: bc.node,
                    dof: bc.dof,
                });
            }
        }
        for l in &loads {
            mesh.dof_index(l.node, l.dof)?;
        }
        Ok(Self {
            mesh,
            bcs,
            loads,
            distributed: Vec::new(),
        })
    }

    pub fn with_distributed_loads(
        mut self,
        distributed: Vec<DistributedLoad>,
    ) -> Result<Self, MailleError> {
        if let Some(d) = distributed
            .iter()
            .find(|d| d.element >= self.mesh.elements.len())
        {
            return Err(MailleError::UnknownElement { element: d.element });
        }
        self.distributed = distributed;
        Ok(self)
    }

    /// Statique linéaire : K_ff u_f = F_f − K_fc u_c, puis réactions R = K u − F sur les DOF contraints.
    /// K n'est jamais formée en entier : K_ff est assemblée en ligne de ciel, les réactions sont
    /// recalculées Element par Element.
    pub fn solve(&self) -> Result<Field, MailleError> {
        let mesh = &self.mesh;
        let n = mesh.nodes.len() * mesh.dofs_per_node();

        let mut f = vec![0.0; n];
        for l in &self.loads {
            f[mesh.dof_index(l.node, l.dof)?] += l.value;
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
            for (g, v) in element_dofs(e).into_iter().zip(equivalent) {
                f[g] += v;
            }
        }
        let mut u = vec![0.0; n];
        let mut fixed = vec![false; n];
        for bc in &self.bcs {
            let i = mesh.dof_index(bc.node, bc.dof)?;
            u[i] = bc.value;
            fixed[i] = true;
        }

        let uf = solve_free(mesh, &fixed, &u, &f)?;
        for (i, v) in (0..n).filter(|&i| !fixed[i]).zip(uf) {
            u[i] = v;
        }

        let mut reactions = vec![0.0; n];
        for i in 0..mesh.elements.len() {
            let dofs = element_dofs(&mesh.elements[i]);
            if !dofs.iter().any(|&g| fixed[g]) {
                continue;
            }
            let ke = element_stiffness(mesh, i);
            for (row, &g) in ke.iter().zip(&dofs).filter(|(_, &g)| fixed[g]) {
                reactions[g] += row.iter().zip(&dofs).map(|(k, &h)| k * u[h]).sum::<f64>();
            }
        }
        for i in (0..n).filter(|&i| fixed[i]) {
            reactions[i] -= f[i];
        }
        let stresses = (0..mesh.elements.len())
            .map(|i| element_stresses(mesh, i, &u))
            .collect();
        Ok(Field {
            displacements: u,
            reactions,
            stresses,
        })
    }
}

/// Résout K_ff u_f = F_f − K_fc u_c. Les DOF libres sont renumérotés par Cuthill-McKee inverse
/// (graphe : deux DOF sont voisins s'ils partagent un Element) avant l'assemblage en ligne de ciel.
/// Renvoie u_f dans l'ordre croissant des DOF libres.
fn solve_free(mesh: &Mesh, fixed: &[bool], u: &[f64], f: &[f64]) -> Result<Vec<f64>, MailleError> {
    // free[g] : indice du DOF global g parmi les DOF libres
    let mut free = vec![None; fixed.len()];
    let mut count = 0;
    for (g, _) in fixed.iter().enumerate().filter(|(_, &c)| !c) {
        free[g] = Some(count);
        count += 1;
    }
    if count == 0 {
        return Ok(Vec::new());
    }
    let element_free: Vec<Vec<usize>> = mesh
        .elements
        .iter()
        .map(|e| element_dofs(e).iter().filter_map(|&g| free[g]).collect())
        .collect();

    let mut adjacency = vec![Vec::new(); count];
    for dofs in &element_free {
        for &a in dofs {
            adjacency[a].extend(dofs.iter().filter(|&&b| b != a));
        }
    }
    for neighbours in &mut adjacency {
        neighbours.sort_unstable();
        neighbours.dedup();
    }
    // rank[a] : position du DOF libre a dans la numérotation du solveur
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
    let mut rhs = vec![0.0; count];
    for (g, a) in free.iter().enumerate() {
        if let Some(a) = a {
            rhs[rank[*a]] += f[g];
        }
    }
    for i in 0..mesh.elements.len() {
        let ke = element_stiffness(mesh, i);
        let dofs = element_dofs(&mesh.elements[i]);
        for (row, &gi) in ke.iter().zip(&dofs) {
            let Some(a) = free[gi] else { continue };
            for (k, &gj) in row.iter().zip(&dofs) {
                match free[gj] {
                    // K symétrique : chaque paire (i, j) est ajoutée une fois, en dessous de la diagonale
                    Some(b) if rank[b] <= rank[a] => kff.add(rank[a], rank[b], *k),
                    Some(_) => {}
                    None => rhs[rank[a]] -= k * u[gj],
                }
            }
        }
    }

    let x = kff.cholesky_solve(&rhs)?;
    Ok(rank.iter().map(|&r| x[r]).collect())
}
