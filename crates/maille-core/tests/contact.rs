//! Contact sans frottement contre une surface rigide, par pénalité : patch test (bloc écrasé sur un
//! plan), décharge, contact de Hertz en 2D (cylindre rigide dans un bloc élastique).

mod common;

use common::{box_grid, rect_grid, to_quadratic, to_quadratic_tet};
use maille_core::contact::{Contact, RigidSurface};
use maille_core::material::Material;
use maille_core::mesh::{Dof, Element, Mesh, PlaneState, Tet10, Tet4, Tri3, Tri6};
use maille_core::mesher::{triangulate, Curve, Geometry, MeshOptions};
use maille_core::model::{Amplitude, Bc, Kinematics, Model, Step};
use maille_core::results::Increment;
use maille_core::MailleError;

const E: f64 = 1000.0;
const NU: f64 = 0.3;
const W: f64 = 2.0;
const H: f64 = 1.5;

fn material() -> Material {
    Material::with_poisson(E, NU).unwrap()
}

fn bc(node: usize, dof: Dof, value: f64) -> Bc {
    Bc { node, dof, value }
}

fn near(a: f64, b: f64) -> bool {
    (a - b).abs() < 1e-9
}

fn plane_meshes() -> Vec<(&'static str, Mesh)> {
    let (nodes, tris) = rect_grid(4, 3, W, H);
    let q = to_quadratic(&nodes, &tris);
    let strain = PlaneState::Strain;
    let tri3 = tris
        .iter()
        .map(|&n| Element::Tri3(Tri3::new(n, material(), 1.0, strain).unwrap()))
        .collect();
    let tri6 = q
        .elements
        .iter()
        .map(|&n| Element::Tri6(Tri6::new(n, material(), 1.0, strain).unwrap()))
        .collect();
    vec![
        ("Tri3", Mesh::new(nodes, tri3).unwrap()),
        ("Tri6", Mesh::new(q.nodes, tri6).unwrap()),
    ]
}

fn solid_meshes() -> Vec<(&'static str, Mesh)> {
    let (nodes, tets) = box_grid([2, 2, 2], [W, W, H]);
    let q = to_quadratic_tet(&nodes, &tets);
    let tet4 = tets
        .iter()
        .map(|&n| Element::Tet4(Tet4::new(n, material())))
        .collect();
    let tet10 = q
        .elements
        .iter()
        .map(|&n| Element::Tet10(Tet10::new(n, material())))
        .collect();
    vec![
        ("Tet4", Mesh::new(nodes, tet4).unwrap()),
        ("Tet10", Mesh::new(q.nodes, tet10).unwrap()),
    ]
}

/// Bloc posé sur le plan z = 0 (y = 0 en 2D), dont la face du dessus descend de δ × λ(t). Le plan
/// rigide, les Nodes du dessous en contact ; un seul Node est bloqué latéralement, le reste glisse.
/// Renvoie (BC, contact, hauteur H, ensemble des Nodes du dessous).
fn block_on_plane(mesh: &Mesh, delta: f64) -> (Vec<Bc>, Contact, Vec<usize>) {
    let solid = mesh.dofs_per_node() == 3;
    let up = if solid { Dof::Uz } else { Dof::Uy };
    let coord = |n: &maille_core::mesh::Node| if solid { n.z } else { n.y };
    let mut bcs = Vec::new();
    let mut bottom = Vec::new();
    for (i, n) in mesh.nodes.iter().enumerate() {
        if near(coord(n), H) {
            bcs.push(bc(i, up, -delta));
        }
        if near(coord(n), 0.0) {
            bottom.push(i);
        }
    }
    // appuis lâches : pas de mouvement de corps rigide dans le plan du dessous
    let pin = |x: f64, y: f64| {
        bottom
            .iter()
            .copied()
            .find(|&i| near(mesh.nodes[i].x, x) && near(mesh.nodes[i].y, y))
            .unwrap()
    };
    let a = pin(0.0, 0.0);
    bcs.push(bc(a, Dof::Ux, 0.0));
    if solid {
        bcs.push(bc(a, Dof::Uy, 0.0));
        bcs.push(bc(pin(W, 0.0), Dof::Uy, 0.0));
    }
    let normal = if solid {
        [0.0, 0.0, 1.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let contact = Contact {
        nodes: bottom.clone(),
        surface: RigidSurface::Plane {
            point: [0.0; 3],
            normal,
        },
        motion: [0.0; 3],
    };
    (bcs, contact, bottom)
}

fn solve(model: &Model, time: f64, increments: usize, amplitude: &Amplitude) -> Vec<Increment> {
    let step = Step::quasi_static(time, increments).unwrap();
    model.solve_quasi_static(&step, amplitude).unwrap()
}

#[test]
fn a_block_crushed_on_a_rigid_plane_is_in_uniform_compression() {
    let delta = 0.03;
    let meshes: Vec<_> = plane_meshes().into_iter().chain(solid_meshes()).collect();
    for (name, mesh) in meshes {
        let solid = mesh.dofs_per_node() == 3;
        let (bcs, contact, bottom) = block_on_plane(&mesh, delta);
        let model = Model::new(mesh.clone(), bcs, Vec::new())
            .unwrap()
            .with_contacts(vec![contact])
            .unwrap();
        let history = solve(&model, 1.0, 3, &Amplitude::ramp(1.0).unwrap());
        let field = &history[3].field;
        let dim = mesh.dofs_per_node();
        // pénétration : la même partout, négligeable devant l'enfoncement
        let up = dim - 1;
        let penetration = -field.displacements[dim * bottom[0] + up];
        for &n in &bottom {
            let u = field.displacements[dim * n + up];
            assert!((u + penetration).abs() < 1e-9 * delta, "{name} : u = {u}");
        }
        assert!(
            penetration > 0.0 && penetration < 1e-3 * delta,
            "{name} : pénétration {penetration} pour un enfoncement {delta}"
        );
        // contrainte uniforme et égale à celle de la déformation réelle (δ − pénétration)
        // côtés libres : σ = E ε (bloc 3D) ou E ε / (1 − ν²) (déformations planes)
        let modulus = if solid { E } else { E / (1.0 - NU * NU) };
        let sigma = -modulus * (delta - penetration) / H;
        for s in field.stresses.iter().flatten() {
            let vertical = if solid { s.sigma_zz } else { s.sigma_yy };
            assert!(
                (vertical - sigma).abs() < 1e-8 * sigma.abs(),
                "{name} : σ = {vertical} au lieu de {sigma}"
            );
        }
        // force sur le plan (rigide, sous le bloc : elle pousse vers le bas) = σ × aire du dessous
        let area = if solid { W * W } else { W };
        let total = field.rigid_forces[0][up];
        assert!(
            (total - sigma * area).abs() < 1e-8 * (sigma * area).abs(),
            "{name} : F = {total} au lieu de {}",
            sigma * area
        );
        // la pression de contact est la contrainte
        for &n in &bottom {
            assert!(
                (field.contact_pressures[n] + sigma).abs() < 1e-6 * sigma.abs(),
                "{name} : p = {}",
                field.contact_pressures[n]
            );
        }
    }
}

#[test]
fn a_moving_plane_does_the_same_as_a_moving_top() {
    // le dessus est fixe (uy = 0) et le plan monte de δ : même état qu'un dessus qui descend de δ
    let delta = 0.03;
    for (name, mesh) in plane_meshes() {
        let (mut bcs, mut contact, _) = block_on_plane(&mesh, delta);
        for b in bcs.iter_mut().filter(|b| b.dof == Dof::Uy) {
            b.value = 0.0;
        }
        contact.motion = [0.0, delta, 0.0];
        let moving = Model::new(mesh.clone(), bcs, Vec::new())
            .unwrap()
            .with_contacts(vec![contact])
            .unwrap();
        let (bcs, contact, _) = block_on_plane(&mesh, delta);
        let fixed = Model::new(mesh, bcs, Vec::new())
            .unwrap()
            .with_contacts(vec![contact])
            .unwrap();
        let ramp = Amplitude::ramp(1.0).unwrap();
        let a = solve(&moving, 1.0, 3, &ramp);
        let b = solve(&fixed, 1.0, 3, &ramp);
        // les déformations sont les mêmes, à la pénétration près : contraintes égales à 1e-3
        let (sa, sb) = (a[3].field.stresses[0][0], b[3].field.stresses[0][0]);
        assert!(
            (sa.sigma_yy - sb.sigma_yy).abs() < 2e-3 * sb.sigma_yy.abs(),
            "{name} : {} et {}",
            sa.sigma_yy,
            sb.sigma_yy
        );
    }
}

#[test]
fn unloading_leaves_no_contact_and_the_block_returns_to_rest() {
    let delta = 0.03;
    for (name, mesh) in plane_meshes().into_iter().chain(solid_meshes()) {
        let (bcs, contact, _) = block_on_plane(&mesh, delta);
        let model = Model::new(mesh, bcs, Vec::new())
            .unwrap()
            .with_contacts(vec![contact])
            .unwrap();
        let amplitude = Amplitude::new(vec![(0.0, 0.0), (1.0, 1.0), (2.0, 0.0)]).unwrap();
        let history = solve(&model, 2.0, 8, &amplitude);
        let loaded = &history[4].field;
        assert!(
            loaded.rigid_forces[0].iter().any(|f| f.abs() > 1e-3),
            "{name}"
        );
        let end = &history[8].field;
        for f in end.contact_forces.iter().chain(&end.contact_pressures) {
            assert!(f.abs() < 1e-10 * E, "{name} : contact résiduel {f}");
        }
        assert!(
            end.rigid_forces[0].iter().all(|f| f.abs() < 1e-10 * E),
            "{name}"
        );
        for u in &end.displacements {
            assert!(u.abs() < 1e-10, "{name} : u = {u}");
        }
    }
}

#[test]
fn contact_works_in_large_strain_too() {
    // même bloc, grandes déformations : Newton converge et la force reste du même ordre qu'en petites
    let delta = 0.03;
    for (name, mesh) in plane_meshes() {
        let (bcs, contact, _) = block_on_plane(&mesh, delta);
        let small = Model::new(mesh.clone(), bcs.clone(), Vec::new())
            .unwrap()
            .with_contacts(vec![contact.clone()])
            .unwrap();
        let large = small.clone().with_kinematics(Kinematics::Large).unwrap();
        let ramp = Amplitude::ramp(1.0).unwrap();
        let (a, b) = (solve(&small, 1.0, 3, &ramp), solve(&large, 1.0, 3, &ramp));
        let (fa, fb) = (a[3].field.rigid_forces[0][1], b[3].field.rigid_forces[0][1]);
        assert!(
            (fa - fb).abs() < 5e-2 * fa.abs(),
            "{name} : {fa} (petites) et {fb} (grandes)"
        );
    }
}

#[test]
fn contact_inputs_are_typed_errors() {
    let (_, mesh) = plane_meshes().remove(0);
    let (bcs, contact, bottom) = block_on_plane(&mesh, 0.01);
    let model = Model::new(mesh.clone(), bcs, Vec::new()).unwrap();
    let with = |c: Contact| model.clone().with_contacts(vec![c]);

    assert_eq!(
        with(Contact {
            nodes: vec![10_000],
            ..contact.clone()
        }),
        Err(MailleError::UnknownNode { node: 10_000 })
    );
    assert!(matches!(
        with(Contact {
            surface: RigidSurface::Plane {
                point: [0.0; 3],
                normal: [0.0; 3]
            },
            ..contact.clone()
        }),
        Err(MailleError::InvalidRigidSurface { .. })
    ));
    assert!(matches!(
        with(Contact {
            surface: RigidSurface::Cylinder {
                point: [0.0; 3],
                axis: [0.0, 0.0, 1.0],
                radius: -1.0
            },
            ..contact.clone()
        }),
        Err(MailleError::InvalidRigidSurface { .. })
    ));
    assert!(matches!(
        with(Contact {
            motion: [f64::NAN, 0.0, 0.0],
            ..contact.clone()
        }),
        Err(MailleError::NonFinite { .. })
    ));
    // un seul Node : aucune face de bord ne le porte
    assert_eq!(
        with(Contact {
            nodes: vec![bottom[0]],
            ..contact.clone()
        }),
        Err(MailleError::ContactNodeWithoutFace { node: bottom[0] })
    );
    // solve linéaire, analyse thermique
    let contacted = with(contact.clone()).unwrap();
    assert_eq!(contacted.solve(), Err(MailleError::ContactNeedsQuasiStatic));
    let thermal = Model::with_analysis(
        mesh,
        maille_core::model::Analysis::Thermal,
        Vec::new(),
        Vec::new(),
    );
    assert!(thermal.is_err() || thermal.unwrap().with_contacts(vec![contact]).is_err());
}

/// Contact de Hertz en 2D : un cylindre rigide de rayon R enfoncé dans un bloc élastique en déformations
/// planes, demi-modèle symétrique en x = 0 (bloc de 10 de côté, 20 fois la demi-largeur de contact).
/// Avec E* = E / (1 − ν²) et P la force par unité d'épaisseur : a = √(4 P R / (π E*)),
/// p0 = 2 P / (π a), p(x) = p0 √(1 − x²/a²). Le déplacement du cylindre n'est pas comparé : en 2D
/// il dépend du bloc (le déplacement d'un demi-plan diverge en log).
mod hertz {
    use super::*;

    pub const R: f64 = 10.0;
    pub const SIZE: f64 = 10.0;

    pub struct Result {
        /// force sur le cylindre, modèle entier
        pub force: f64,
        /// demi-largeur de contact lue sur les Nodes en contact
        pub half_width: f64,
        /// plus grande pression de contact aux Nodes
        pub peak: f64,
        pub penetration: f64,
    }

    pub fn solve(fine: f64, delta: f64, increments: usize) -> Result {
        let geometry = Geometry {
            outer: vec![
                Curve::line([0.0, 0.0], fine, "contact"),
                Curve::line([1.5, 0.0], 8.0 * fine, "surface"),
                Curve::line([SIZE, 0.0], 2.0, "droite"),
                Curve::line([SIZE, -SIZE], 2.0, "bas"),
                Curve::line([0.0, -SIZE], 1.0, "sym"),
            ],
            holes: vec![],
        };
        let t = triangulate(&geometry, &MeshOptions::default())
            .unwrap()
            .to_quadratic();
        let node_set = |name: &str| t.boundary(name).unwrap().nodes.clone();
        let (top1, top2) = (node_set("contact"), node_set("surface"));
        let (sym, bottom, right) = (node_set("sym"), node_set("bas"), node_set("droite"));
        let elements = t
            .triangles
            .iter()
            .map(|&n| Element::Tri6(Tri6::new(n, material(), 1.0, PlaneState::Strain).unwrap()))
            .collect();
        let mesh = Mesh::new(t.nodes, elements).unwrap();

        let mut bcs: Vec<Bc> = sym.iter().map(|&n| bc(n, Dof::Ux, 0.0)).collect();
        for &n in bottom.iter().chain(&right) {
            bcs.push(bc(n, Dof::Ux, 0.0));
            bcs.push(bc(n, Dof::Uy, 0.0));
        }
        // le Node commun du bas et du côté a deux BC identiques : on garde la première
        let mut seen = std::collections::HashSet::new();
        bcs.retain(|b| seen.insert((b.node, b.dof as u8)));

        let mut slave: Vec<usize> = top1.iter().chain(&top2).copied().collect();
        slave.sort_unstable();
        slave.dedup();
        let contact = Contact {
            nodes: slave.clone(),
            surface: RigidSurface::Cylinder {
                point: [0.0, R, 0.0],
                axis: [0.0, 0.0, 1.0],
                radius: R,
            },
            motion: [0.0, -delta, 0.0],
        };
        let model = Model::new(mesh.clone(), bcs, Vec::new())
            .unwrap()
            .with_contacts(vec![contact])
            .unwrap();
        let history = super::solve(&model, 1.0, increments, &Amplitude::ramp(1.0).unwrap());
        let field = &history[increments].field;

        let mut profile: Vec<(f64, f64)> = slave
            .iter()
            .map(|&n| (mesh.nodes[n].x, field.contact_pressures[n]))
            .filter(|(_, p)| *p > 0.0)
            .collect();
        profile.sort_by(|a, b| a.0.total_cmp(&b.0));
        let last = profile.last().unwrap().0;
        // le Node suivant, hors contact : la largeur est entre les deux
        let next = slave
            .iter()
            .map(|&n| mesh.nodes[n].x)
            .filter(|&x| x > last + 1e-9)
            .fold(f64::INFINITY, f64::min);
        Result {
            force: 2.0 * field.rigid_forces[0][1],
            half_width: 0.5 * (last + next),
            peak: profile.iter().fold(0.0, |m, p| f64::max(m, p.1)),
            // g = distance au centre − R, le Node monte de uy et le centre descend de δ
            penetration: slave
                .iter()
                .map(|&n| {
                    let x = mesh.nodes[n].x;
                    let y = field.displacements[2 * n + 1];
                    R - x.hypot(y - (R - delta))
                })
                .fold(0.0, f64::max),
        }
    }
}

#[test]
fn hertz_2d_matches_the_half_width_and_the_peak_pressure() {
    let delta = 0.08;
    let r = hertz::solve(0.04, delta, 8);
    let e_star = E / (1.0 - NU * NU);
    let a = (4.0 * r.force * hertz::R / (std::f64::consts::PI * e_star)).sqrt();
    let p0 = 2.0 * r.force / (std::f64::consts::PI * a);
    assert!(
        (r.half_width / a - 1.0).abs() < 5e-2,
        "a = {}, Hertz {a}",
        r.half_width
    );
    assert!(
        (r.peak / p0 - 1.0).abs() < 5e-2,
        "p0 = {}, Hertz {p0}",
        r.peak
    );
    assert!(
        r.penetration < 1e-3 * delta,
        "pénétration {}",
        r.penetration
    );
}
