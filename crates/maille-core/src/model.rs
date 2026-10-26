use std::ops::Deref;

use crate::assembly::{
    element_capacity, element_conduction, element_dofs, element_gauss_count,
    element_heat_flux_operator, element_heat_source, element_response, element_stiffness,
    element_stresses, element_thermal_load, element_thermal_strains, face_convection,
    face_pressure, FaceConvection, FluxOperator,
};
use crate::element::{bar_equivalent_load, beam_equivalent_load, ElementResponse};
use crate::error::MailleError;
use crate::material::positive;
use crate::mesh::{Dof, Element, Mesh, PlaneState};
use crate::plasticity::PlasticState;
use crate::results::{Field, Increment};
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

/// Pression uniforme sur une face d'un Element (statique). Numérotation des faces : voir
/// `Element::face_nodes`. Positive, elle comprime (elle agit contre la normale sortante, comme
/// `*DLOAD P` d'Abaqus) ; négative, elle tire.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Pressure {
    pub element: usize,
    pub face: usize,
    /// Une contrainte (ex. MPa), finie.
    pub value: f64,
}

/// Physique résolue sur le Mesh.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Analysis {
    /// Statique : déplacements, contraintes ; linéaire avec `Model::solve`, élasto-plastique (et p)
    /// avec `Model::solve_quasi_static`.
    Static,
    /// Thermique, un DOF `T` par Node : stationnaire −∇·(k ∇T) = Q avec `Model::solve`, transitoire
    /// ρc ∂T/∂t − ∇·(k ∇T) = Q avec `Model::solve_step`.
    Thermal,
}

/// Step : durée découpée en incréments égaux Δt. Deux usages :
/// - transitoire (`Step::transient`, Model `Thermal`, voir `Model::solve_step`) : θ-méthode
///   C (Tₙ₊₁ − Tₙ)/Δt + K (θ Tₙ₊₁ + (1 − θ) Tₙ) = F ; les BC, sources et convections du Model restent
///   constantes pendant la Step ;
/// - quasi-statique (`Step::quasi_static`, Model `Static`, voir `Model::solve_quasi_static`) :
///   charges et valeurs imposées des BC suivent une `Amplitude` ; le temps ne sert qu'à la lire.
///
/// La Step ne porte que le temps : l'Analysis reste sur le Model, y compris en couplage séquentiel
/// (voir `with_temperatures`).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Step {
    time_period: f64,
    increments: usize,
    theta: f64,
}

impl Step {
    /// Durée `time_period` (finie, strictement positive, dans l'unité de temps du modèle) en
    /// `increments` pas égaux, en Euler implicite (θ = 1).
    pub fn transient(time_period: f64, increments: usize) -> Result<Self, MailleError> {
        let time_period = positive("durée de la Step", time_period)?;
        if increments == 0 {
            return Err(MailleError::NoIncrements);
        }
        Ok(Self {
            time_period,
            increments,
            theta: 1.0,
        })
    }

    /// Durée `time_period` (finie, strictement positive) en `increments` pas égaux, pour un calcul
    /// élasto-plastique : chaque incrément est un équilibre résolu par Newton-Raphson. Plus
    /// d'incréments aident Newton à converger quand la plasticité se développe vite. C'est la même
    /// Step que `transient`, nommée pour son usage : θ (`with_theta`) n'y a aucun effet.
    pub fn quasi_static(time_period: f64, increments: usize) -> Result<Self, MailleError> {
        Self::transient(time_period, increments)
    }

    /// θ entre ½ et 1, stable quel que soit Δt. 1 : Euler implicite, d'ordre 1, qui amortit tout ;
    /// à préférer après un départ discontinu (une trempe). ½ : Crank-Nicolson, d'ordre 2, mais qui
    /// oscille après un départ discontinu si Δt est grand.
    pub fn with_theta(self, theta: f64) -> Result<Self, MailleError> {
        if !(0.5..=1.0).contains(&theta) {
            return Err(MailleError::ThetaOutOfRange { theta });
        }
        Ok(Self { theta, ..self })
    }

    pub fn time_period(&self) -> f64 {
        self.time_period
    }

    pub fn increments(&self) -> usize {
        self.increments
    }

    pub fn theta(&self) -> f64 {
        self.theta
    }

    /// Pas de temps Δt.
    pub fn time_step(&self) -> f64 {
        self.time_period / self.increments as f64
    }
}

/// Amplitude λ(t) d'une Step quasi-statique : les charges (Load, DistributedLoad) et les valeurs
/// des BC du Model sont multipliées par λ(t). Linéaire par morceaux entre des points (t, λ),
/// constante après le dernier (comme *AMPLITUDE dans Abaqus). Charger puis décharger :
/// `Amplitude::new(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)])`.
#[derive(Debug, Clone, PartialEq)]
pub struct Amplitude {
    points: Vec<(f64, f64)>,
}

impl Amplitude {
    /// Le premier point doit être (0, 0) : le calcul part du Mesh au repos. Les temps sont
    /// strictement croissants et les valeurs finies.
    pub fn new(points: Vec<(f64, f64)>) -> Result<Self, MailleError> {
        if points.is_empty() {
            return Err(MailleError::InvalidAmplitude { point: 0 });
        }
        for (k, &(t, value)) in points.iter().enumerate() {
            let valid = t.is_finite()
                && value.is_finite()
                && match k {
                    // (0, 0) à l'arrondi près
                    0 => t.abs() <= f64::EPSILON && value.abs() <= f64::EPSILON,
                    _ => t > points[k - 1].0,
                };
            if !valid {
                return Err(MailleError::InvalidAmplitude { point: k });
            }
        }
        Ok(Self { points })
    }

    /// Rampe de 0 à 1 sur `time_period`, puis 1.
    pub fn ramp(time_period: f64) -> Result<Self, MailleError> {
        Self::new(vec![
            (0.0, 0.0),
            (positive("durée de la rampe", time_period)?, 1.0),
        ])
    }

    /// λ(t).
    pub fn at(&self, t: f64) -> f64 {
        let after = self.points.partition_point(|&(ti, _)| ti <= t);
        match (
            self.points.get(after.wrapping_sub(1)),
            self.points.get(after),
        ) {
            (Some(&(t0, v0)), Some(&(t1, v1))) => v0 + (v1 - v0) * (t - t0) / (t1 - t0),
            (Some(&(_, v)), None) | (None, Some(&(_, v))) => v,
            (None, None) => 0.0,
        }
    }
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
    pressures: Vec<Pressure>,
    /// Champ de températures aux Nodes et T_ref (statique), voir `with_temperatures`.
    temperatures: Option<(Vec<f64>, f64)>,
}

impl Model {
    /// Analyse statique, voir `with_analysis`.
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
                Analysis::Static if m.plasticity.is_some() && !elastoplastic(e) => {
                    return Err(MailleError::UnsupportedPlasticElement { element: i })
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
            pressures: Vec::new(),
            temperatures: None,
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

    /// Pressions sur des faces, en statique seulement : triangles (Tri3, Tri6) et tétraèdres (Tet4,
    /// Tet10). Charge morte en petites déformations : elle ne suit pas la déformée. En quasi-statique,
    /// elle suit l'Amplitude comme les Load. En 2D, l'épaisseur de l'Element multiplie l'intégrale.
    pub fn with_pressures(mut self, pressures: Vec<Pressure>) -> Result<Self, MailleError> {
        if self.analysis != Analysis::Static {
            return Err(MailleError::PressureStaticOnly);
        }
        for p in &pressures {
            let e = self.element(p.element)?;
            if matches!(e, Element::Bar(_) | Element::Beam(_)) {
                return Err(MailleError::UnsupportedPressureElement { element: p.element });
            }
            if e.face_positions(p.face).is_none() {
                return Err(MailleError::UnknownFace {
                    element: p.element,
                    face: p.face,
                });
            }
            if !p.value.is_finite() {
                return Err(MailleError::NonFinite {
                    name: "p (pression)",
                    value: p.value,
                });
            }
        }
        self.pressures = pressures;
        Ok(self)
    }

    /// Champ de températures qui dilate le Mesh, en statique seulement (couplage thermo-mécanique
    /// séquentiel) : une valeur par Node, par exemple les `temperatures` d'un Field thermique calculé
    /// sur le même Mesh. Ce n'est pas une BC : aucun DOF `T` n'est résolu en statique.
    /// Chaque Element se dilate de α (T − `reference`), `reference` étant la température sans
    /// contrainte ; les contraintes du Field sont alors σ = D (ε − ε_th). Bar, triangles et
    /// tétraèdres, dont le Material doit avoir α.
    pub fn with_temperatures(
        mut self,
        temperatures: Vec<f64>,
        reference: f64,
    ) -> Result<Self, MailleError> {
        if self.analysis != Analysis::Static {
            return Err(MailleError::StaticOnly);
        }
        check_temperatures(&self.mesh, &temperatures)?;
        if !reference.is_finite() {
            return Err(MailleError::NonFinite {
                name: "T_ref (température sans contrainte)",
                value: reference,
            });
        }
        for (i, e) in self.mesh.elements.iter().enumerate() {
            if e.material().plasticity.is_some() {
                return Err(MailleError::QuasiStaticWithTemperatures);
            }
            if matches!(e, Element::Beam(_)) {
                return Err(MailleError::UnsupportedThermalStrain { element: i });
            }
            if e.material().expansion.is_none() {
                return Err(MailleError::MissingExpansion { element: i });
            }
        }
        self.temperatures = Some((temperatures, reference));
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
        if let Some(i) = mesh
            .elements
            .iter()
            .position(|e| e.material().plasticity.is_some())
        {
            return Err(MailleError::PlasticityNeedsQuasiStatic { element: i });
        }
        let mut a = self.assemble()?;
        // dilatation : α (T − T_ref) aux Nodes de chaque Element (vide sans champ de températures),
        // et ses forces nodales dans F
        let mut strains = vec![Vec::new(); mesh.elements.len()];
        if let Some((t, reference)) = &self.temperatures {
            for (i, strain_e) in strains.iter_mut().enumerate() {
                *strain_e = element_thermal_strains(mesh, i, t, *reference);
                for (&g, v) in a.dofs[i]
                    .iter()
                    .zip(element_thermal_load(mesh, i, strain_e))
                {
                    a.f[g] += v;
                }
            }
        }
        let matrix = |i: usize| match self.analysis {
            Analysis::Static => element_stiffness(mesh, i),
            Analysis::Thermal => thermal_matrix(mesh, &a.faces, i),
        };
        let mut u = a.imposed.clone();
        FreeSystem::new(&a.dofs, &matrix, &a.fixed, &u)?.solve(&a.f, &mut u);

        let mut reactions = product(&a.dofs, &matrix, &u, &a.fixed);
        for (g, r) in reactions
            .iter_mut()
            .enumerate()
            .filter(|(g, _)| a.fixed[*g])
        {
            *r -= a.f[g];
        }
        Ok(match self.analysis {
            Analysis::Static => Field {
                stresses: (0..mesh.elements.len())
                    .map(|i| element_stresses(mesh, i, &u, &strains[i]))
                    .collect(),
                equivalent_plastic_strains: Vec::new(),
                displacements: u,
                reactions,
                temperatures: Vec::new(),
                heat_fluxes: Vec::new(),
            },
            Analysis::Thermal => thermal_field(&flux_operators(mesh), u, reactions),
        })
    }

    /// Thermique transitoire : intègre la Step depuis `initial`, une température par Node (par
    /// exemple `vec![T0; nombre de Nodes]`, ou les `temperatures` d'un Field stationnaire). Les BC
    /// s'appliquent dès t = 0 : elles remplacent la condition initiale à leurs Nodes. Renvoie
    /// l'incrément t = 0, puis un par pas.
    ///
    /// À chaque pas, (C/Δt + θK) Tₙ₊₁ = (C/Δt − (1 − θ)K) Tₙ + F, K comprenant la convection. Δt est
    /// fixe : la matrice de gauche est factorisée une seule fois. Réactions sur les Nodes imposés :
    /// R = C (Tₙ₊₁ − Tₙ)/Δt + K (θ Tₙ₊₁ + (1 − θ) Tₙ) − F, puissance moyenne fournie pendant le pas
    /// (Σ R Δt est exactement la chaleur apportée par les BC).
    pub fn solve_step(&self, step: &Step, initial: &[f64]) -> Result<Vec<Increment>, MailleError> {
        self.require_thermal()?;
        let mesh = &self.mesh;
        let n = mesh.nodes.len();
        check_temperatures(mesh, initial)?;
        for (i, e) in mesh.elements.iter().enumerate() {
            if e.material().density.is_none() {
                return Err(MailleError::MissingDensity { element: i });
            }
            if e.material().specific_heat.is_none() {
                return Err(MailleError::MissingSpecificHeat { element: i });
            }
        }

        let a = self.assemble()?;
        let (dt, theta) = (step.time_step(), step.theta());
        // matrices de chaque Element : à gauche C/Δt + θK, à droite C/Δt − (1 − θ)K
        let (left, right): (Vec<_>, Vec<_>) = (0..mesh.elements.len())
            .map(|i| {
                let (c, k) = (element_capacity(mesh, i), thermal_matrix(mesh, &a.faces, i));
                let combine = |s: f64| -> Vec<Vec<f64>> {
                    c.iter()
                        .zip(&k)
                        .map(|(cr, kr)| cr.iter().zip(kr).map(|(c, k)| c / dt + s * k).collect())
                        .collect()
                };
                (combine(theta), combine(theta - 1.0))
            })
            .unzip();
        let (left, right) = (
            |i: usize| left[i].as_slice(),
            |i: usize| right[i].as_slice(),
        );

        let mut t = initial.to_vec();
        for (g, _) in a.fixed.iter().enumerate().filter(|(_, &c)| c) {
            t[g] = a.imposed[g];
        }
        let system = FreeSystem::new(&a.dofs, &left, &a.fixed, &t)?;
        let all = vec![true; n];
        let fluxes = flux_operators(mesh);
        let mut history = vec![Increment {
            time: 0.0,
            field: thermal_field(&fluxes, t.clone(), vec![0.0; n]),
            residuals: Vec::new(),
        }];
        for s in 1..=step.increments() {
            let carried = product(&a.dofs, &right, &t, &all);
            let rhs: Vec<f64> = carried.iter().zip(&a.f).map(|(c, f)| c + f).collect();
            system.solve(&rhs, &mut t);
            let mut reactions = product(&a.dofs, &left, &t, &a.fixed);
            for (g, r) in reactions
                .iter_mut()
                .enumerate()
                .filter(|(g, _)| a.fixed[*g])
            {
                *r -= rhs[g];
            }
            history.push(Increment {
                time: s as f64 * dt,
                field: thermal_field(&fluxes, t.clone(), reactions),
                residuals: Vec::new(),
            });
        }
        Ok(history)
    }

    /// Calcul quasi-statique, élasto-plastique si des Materials ont une `Plasticity` : les charges et
    /// les valeurs des BC du Model sont multipliées par λ(t) de l'`amplitude`, et chaque incrément de
    /// la Step cherche l'équilibre F_int(u) = λ F par Newton-Raphson, avec la tangente cohérente de la
    /// loi de von Mises (convergence quadratique). L'état plastique de chaque point de Gauss n'est
    /// gardé qu'à la convergence. Renvoie l'incrément t = 0 (Mesh au repos), puis un par pas.
    /// Réactions : F_int − λ F aux DOF imposés.
    ///
    /// Triangles en déformations planes et tétraèdres seulement, sans champ de températures. Un
    /// Material sans `Plasticity` reste élastique : le résultat est alors celui de `solve` mis à
    /// l'échelle par λ, chaque incrément convergeant en une itération.
    pub fn solve_quasi_static(
        &self,
        step: &Step,
        amplitude: &Amplitude,
    ) -> Result<Vec<Increment>, MailleError> {
        if self.analysis != Analysis::Static {
            return Err(MailleError::StaticOnly);
        }
        if self.temperatures.is_some() {
            return Err(MailleError::QuasiStaticWithTemperatures);
        }
        let mesh = &self.mesh;
        if let Some(i) = mesh.elements.iter().position(|e| !elastoplastic(e)) {
            return Err(MailleError::UnsupportedPlasticElement { element: i });
        }
        let a = self.assemble()?;
        let mut committed: Vec<Vec<PlasticState>> = (0..mesh.elements.len())
            .map(|i| vec![PlasticState::default(); element_gauss_count(mesh, i)])
            .collect();
        let mut u = vec![0.0; a.f.len()];
        let at_rest = self.responses(&u, &committed);
        let mut history = vec![Increment {
            time: 0.0,
            field: quasi_static_field(u.clone(), vec![0.0; a.f.len()], at_rest),
            residuals: Vec::new(),
        }];
        // échelle des forces : la plus grande norme de F_int ou de λ F vue depuis le début de la Step
        // (F_int contient les réactions, la référence reste non nulle sous déplacement imposé et
        // quand la décharge ramène les forces à zéro)
        let mut scale = 0.0;
        for increment in 1..=step.increments() {
            let time = increment as f64 * step.time_step();
            let lambda = amplitude.at(time);
            let (responses, reactions, residuals) = self.equilibrium(
                &a,
                &committed,
                &mut u,
                lambda,
                &mut scale,
                (increment, time),
            )?;
            committed = responses.iter().map(|r| r.states.clone()).collect();
            history.push(Increment {
                time,
                field: quasi_static_field(u.clone(), reactions, responses),
                residuals,
            });
        }
        Ok(history)
    }

    /// Newton-Raphson d'un incrément : part de `u` (équilibre de l'incrément précédent, états
    /// `committed`) et le mène à F_int(u) = λ F, valeurs imposées λ × valeur de la BC. Renvoie les
    /// réponses des Elements, les réactions et les résidus relatifs (voir `Increment::residuals`).
    /// `at` : (incrément, temps), pour les erreurs.
    fn equilibrium(
        &self,
        a: &Assembly,
        committed: &[Vec<PlasticState>],
        u: &mut [f64],
        lambda: f64,
        scale: &mut f64,
        (increment, time): (usize, f64),
    ) -> Result<Equilibrium, MailleError> {
        let external: Vec<f64> = a.f.iter().map(|f| lambda * f).collect();
        let mut residuals = Vec::new();
        loop {
            let responses = self.responses(u, committed);
            let internal = internal_force(&responses, &a.dofs, u.len());
            let residual: Vec<f64> = internal.iter().zip(&external).map(|(i, e)| i - e).collect();
            *scale = scale
                .max(norm(internal.iter().copied()))
                .max(norm(external.iter().copied()));
            let free = norm(
                residual
                    .iter()
                    .zip(&a.fixed)
                    .filter(|(_, &c)| !c)
                    .map(|(r, _)| *r),
            );
            let relative = if *scale > 0.0 { free / *scale } else { free };
            if !relative.is_finite() {
                return Err(MailleError::NotConverged { increment, time });
            }
            residuals.push(relative);
            // au moins une itération : la première applique l'incrément des valeurs imposées
            if residuals.len() > 1 && relative <= NEWTON_TOLERANCE {
                let reactions = residual
                    .iter()
                    .zip(&a.fixed)
                    .map(|(r, &c)| if c { *r } else { 0.0 })
                    .collect();
                return Ok((responses, reactions, residuals));
            }
            if residuals.len() > NEWTON_ITERATIONS {
                return Err(MailleError::NotConverged { increment, time });
            }
            // δu impose aux DOF contraints ce qui manque pour atteindre λ × valeur de la BC (non
            // nul à la première itération seulement) : K_ff δu_f = −r_f − K_fc δu_c. Partir de
            // l'état convergé plutôt que de déplacer d'un coup les Nodes imposés évite une
            // plastification artificielle des Elements du bord.
            let mut du: Vec<f64> = a
                .fixed
                .iter()
                .zip(&a.imposed)
                .zip(&*u)
                .map(|((&c, imposed), u)| if c { lambda * imposed - u } else { 0.0 })
                .collect();
            let tangent = |i: usize| responses[i].tangent.as_slice();
            let system = match FreeSystem::new(&a.dofs, &tangent, &a.fixed, &du) {
                Err(MailleError::Singular) if responses.iter().any(|r| r.yielding) => {
                    return Err(MailleError::PlasticCollapse { increment, time })
                }
                other => other?,
            };
            system.solve(&residual.iter().map(|r| -r).collect::<Vec<_>>(), &mut du);
            for (u, du) in u.iter_mut().zip(du) {
                *u += du;
            }
        }
    }

    /// Réponse élasto-plastique de chaque Element aux déplacements `u`, depuis les états `committed`.
    fn responses(&self, u: &[f64], committed: &[Vec<PlasticState>]) -> Vec<ElementResponse> {
        committed
            .iter()
            .enumerate()
            .map(|(i, previous)| element_response(&self.mesh, i, u, previous))
            .collect()
    }

    /// Ce que les résolutions partagent : DOF des Elements, convection, F et valeurs imposées.
    fn assemble(&self) -> Result<Assembly, MailleError> {
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
        for p in &self.pressures {
            if let Some(terms) = face_pressure(mesh, p.element, p.face, p.value) {
                terms.add_load(&dofs[p.element], &mut f);
            }
        }
        for (terms, element_dofs) in faces.iter().zip(&dofs) {
            for face in terms {
                face.add_load(element_dofs, &mut f);
            }
        }
        let mut imposed = vec![0.0; n];
        let mut fixed = vec![false; n];
        for bc in &self.bcs {
            let i = dof_index(mesh, self.analysis, bc.node, bc.dof)?;
            imposed[i] = bc.value;
            fixed[i] = true;
        }
        Ok(Assembly {
            dofs,
            faces,
            f,
            imposed,
            fixed,
        })
    }
}

/// Équilibre d'un incrément : réponses des Elements, réactions, résidus relatifs.
type Equilibrium = (Vec<ElementResponse>, Vec<f64>, Vec<f64>);

/// Tolérance de Newton-Raphson sur le résidu relatif (voir `Increment::residuals`).
const NEWTON_TOLERANCE: f64 = 1e-8;
/// Nombre maximal d'itérations de Newton-Raphson par incrément.
const NEWTON_ITERATIONS: usize = 25;

/// Triangle en déformations planes ou tétraèdre : les Elements du calcul élasto-plastique.
fn elastoplastic(e: &Element) -> bool {
    match e {
        Element::Tri3(t) => t.state == PlaneState::Strain,
        Element::Tri6(t) => t.state == PlaneState::Strain,
        Element::Tet4(_) | Element::Tet10(_) => true,
        Element::Bar(_) | Element::Beam(_) => false,
    }
}

fn norm(v: impl Iterator<Item = f64>) -> f64 {
    v.map(|x| x * x).sum::<f64>().sqrt()
}

/// Forces internes globales Σ_e F_int,e, sur `n` DOF.
fn internal_force(responses: &[ElementResponse], dofs: &[Vec<usize>], n: usize) -> Vec<f64> {
    let mut f = vec![0.0; n];
    for (r, dofs) in responses.iter().zip(dofs) {
        for (&g, v) in dofs.iter().zip(&r.internal_force) {
            f[g] += v;
        }
    }
    f
}

/// Field d'un incrément quasi-statique : contraintes et p aux points de Gauss des `responses`.
fn quasi_static_field(
    displacements: Vec<f64>,
    reactions: Vec<f64>,
    responses: Vec<ElementResponse>,
) -> Field {
    let (stresses, equivalent_plastic_strains) = responses
        .into_iter()
        .map(|r| (r.stresses, r.states.iter().map(|s| s.p).collect()))
        .unzip();
    Field {
        displacements,
        reactions,
        stresses,
        equivalent_plastic_strains,
        temperatures: Vec::new(),
        heat_fluxes: Vec::new(),
    }
}

/// Une température finie par Node.
fn check_temperatures(mesh: &Mesh, t: &[f64]) -> Result<(), MailleError> {
    if t.len() != mesh.nodes.len() {
        return Err(MailleError::TemperatureCount {
            expected: mesh.nodes.len(),
            got: t.len(),
        });
    }
    match t.iter().position(|t| !t.is_finite()) {
        Some(node) => Err(MailleError::NonFiniteTemperature { node }),
        None => Ok(()),
    }
}

/// Données communes aux résolutions (`solve`, `solve_step`, `solve_quasi_static`).
struct Assembly {
    /// DOF globaux de chaque Element, dans l'ordre de sa matrice
    dofs: Vec<Vec<usize>>,
    /// termes de convection de chaque Element
    faces: Vec<Vec<FaceConvection>>,
    f: Vec<f64>,
    /// valeur imposée aux DOF `fixed`, zéro ailleurs
    imposed: Vec<f64>,
    fixed: Vec<bool>,
}

/// Matrice thermique de l'Element `i` : conduction plus convection sur ses faces.
fn thermal_matrix(mesh: &Mesh, faces: &[Vec<FaceConvection>], i: usize) -> Vec<Vec<f64>> {
    let mut k = element_conduction(mesh, i);
    for face in &faces[i] {
        face.add_matrix(&mut k);
    }
    k
}

/// Opérateurs de flux de chaque Element, voir `element_heat_flux_operator`.
fn flux_operators(mesh: &Mesh) -> Vec<FluxOperator> {
    (0..mesh.elements.len())
        .map(|i| element_heat_flux_operator(mesh, i))
        .collect()
}

/// Field thermique : températures, réactions et flux aux points de Gauss (`fluxes` : un opérateur
/// par Element).
fn thermal_field(fluxes: &[FluxOperator], temperatures: Vec<f64>, reactions: Vec<f64>) -> Field {
    Field {
        heat_fluxes: fluxes.iter().map(|flux| flux(&temperatures)).collect(),
        temperatures,
        reactions,
        displacements: Vec::new(),
        stresses: Vec::new(),
        equivalent_plastic_strains: Vec::new(),
    }
}

/// (Σ_e M_e x_e) sur les DOF `rows`, zéro ailleurs ; `matrix(i)` : matrice M_e de l'Element i, DOF
/// `dofs[i]`. Les Elements sans DOF dans `rows` sont sautés.
fn product<M: Deref<Target = [Vec<f64>]>>(
    dofs: &[Vec<usize>],
    matrix: &impl Fn(usize) -> M,
    x: &[f64],
    rows: &[bool],
) -> Vec<f64> {
    let mut y = vec![0.0; rows.len()];
    for (i, dofs) in dofs.iter().enumerate() {
        if !dofs.iter().any(|&g| rows[g]) {
            continue;
        }
        for (row, &g) in matrix(i).iter().zip(dofs).filter(|(_, &g)| rows[g]) {
            y[g] += row.iter().zip(dofs).map(|(k, &h)| k * x[h]).sum::<f64>();
        }
    }
    y
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
    fn new<M: Deref<Target = [Vec<f64>]>>(
        dofs: &[Vec<usize>],
        matrix: &impl Fn(usize) -> M,
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
