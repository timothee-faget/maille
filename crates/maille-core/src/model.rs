use crate::assembly::{assemble, element_dofs, element_stresses, DenseMatrix};
use crate::element::{bar_equivalent_load, beam_equivalent_load};
use crate::error::MailleError;
use crate::mesh::{Dof, Element, Mesh};
use crate::results::Field;
use crate::solve::cholesky_solve;

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
    pub fn solve(&self) -> Result<Field, MailleError> {
        let mesh = &self.mesh;
        let n = mesh.nodes.len() * mesh.dofs_per_node();
        let k = assemble(mesh);

        let mut f = vec![0.0; n];
        for l in &self.loads {
            f[mesh.dof_index(l.node, l.dof)?] += l.value;
        }
        for d in &self.distributed {
            let e = &mesh.elements[d.element];
            let equivalent = match e {
                Element::Bar(_) => bar_equivalent_load(d.value, mesh.length(d.element)).to_vec(),
                Element::Beam(_) => beam_equivalent_load(d.value, mesh.length(d.element)).to_vec(),
                Element::Tri3(_) | Element::Tri6(_) => {
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

        let free: Vec<usize> = (0..n).filter(|&i| !fixed[i]).collect();
        if !free.is_empty() {
            let mut kff = DenseMatrix::zeros(free.len());
            let mut rhs = vec![0.0; free.len()];
            for (a, &i) in free.iter().enumerate() {
                rhs[a] = f[i];
                for j in 0..n {
                    if fixed[j] {
                        rhs[a] -= k.get(i, j) * u[j];
                    }
                }
                for (b, &j) in free.iter().enumerate() {
                    kff.add(a, b, k.get(i, j));
                }
            }
            let uf = cholesky_solve(&kff, &rhs)?;
            for (a, &i) in free.iter().enumerate() {
                u[i] = uf[a];
            }
        }

        let mut reactions = vec![0.0; n];
        for i in (0..n).filter(|&i| fixed[i]) {
            reactions[i] = (0..n).map(|j| k.get(i, j) * u[j]).sum::<f64>() - f[i];
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
