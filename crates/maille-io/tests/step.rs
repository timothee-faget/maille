//! Import STEP : topologie B-rep, surfaces des faces et courbes des arêtes. Pièces écrites à la
//! main (pavé, quart de cylindre à arcs en CIRCLE et en B-spline rationnelle, cylindre à couture
//! comme chez OCCT, tronc de cône en degrés), validées par les longueurs exactes des arêtes, la
//! formule d'Euler, l'usage des arêtes et le volume des tétraèdres ; unités et refus avec message.

use std::collections::HashMap;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

use maille_core::mesher::{tetrahedralize_brep, TetOptions};
use maille_io::step::{read_brep, Brep, EdgeCurve, FaceSurface};

/// Fichier STEP en construction : une instance par ligne, numérotées dans l'ordre.
#[derive(Default)]
struct Step {
    lines: Vec<String>,
}

impl Step {
    fn add(&mut self, entity: &str) -> u64 {
        self.lines.push(entity.to_string());
        self.lines.len() as u64
    }

    fn point(&mut self, p: [f64; 3]) -> u64 {
        self.add(&format!(
            "CARTESIAN_POINT('',({:?},{:?},{:?}))",
            p[0], p[1], p[2]
        ))
    }

    fn direction(&mut self, d: [f64; 3]) -> u64 {
        self.add(&format!("DIRECTION('',({:?},{:?},{:?}))", d[0], d[1], d[2]))
    }

    fn placement(&mut self, origin: [f64; 3]) -> u64 {
        self.frame(origin, [0.0, 0.0, 1.0], [1.0, 0.0, 0.0])
    }

    /// Repère d'origine `origin`, d'axe `z` et de direction de référence `x`.
    fn frame(&mut self, origin: [f64; 3], z: [f64; 3], x: [f64; 3]) -> u64 {
        let (o, z, x) = (self.point(origin), self.direction(z), self.direction(x));
        self.add(&format!("AXIS2_PLACEMENT_3D('',#{o},#{z},#{x})"))
    }

    fn plane(&mut self, origin: [f64; 3], normal: [f64; 3]) -> u64 {
        let x = if normal[0].abs() < 0.9 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 1.0, 0.0]
        };
        let frame = self.frame(origin, normal, x);
        self.add(&format!("PLANE('',#{frame})"))
    }

    /// Cylindre d'axe z.
    fn cylinder(&mut self, radius: f64) -> u64 {
        let frame = self.placement([0.0; 3]);
        self.add(&format!("CYLINDRICAL_SURFACE('',#{frame},{radius:?})"))
    }

    fn vertex(&mut self, p: [f64; 3]) -> u64 {
        let point = self.point(p);
        self.add(&format!("VERTEX_POINT('',#{point})"))
    }

    fn edge(&mut self, v1: u64, v2: u64, curve: u64) -> u64 {
        self.add(&format!("EDGE_CURVE('',#{v1},#{v2},#{curve},.T.)"))
    }

    fn line(&mut self, v1: u64, v2: u64, a: [f64; 3], b: [f64; 3]) -> u64 {
        let p = self.point(a);
        let d = self.direction([b[0] - a[0], b[1] - a[1], b[2] - a[2]]);
        let vector = self.add(&format!("VECTOR('',#{d},1.)"));
        let line = self.add(&format!("LINE('',#{p},#{vector})"));
        self.edge(v1, v2, line)
    }

    /// Face nommée `name` sur la surface `surface` : un contour extérieur d'arêtes (`true` : dans
    /// leur sens).
    fn face(&mut self, name: &str, edges: &[(u64, bool)], surface: u64) -> u64 {
        let oriented: Vec<String> = edges
            .iter()
            .map(|&(e, forward)| {
                let o = self.add(&format!(
                    "ORIENTED_EDGE('',*,*,#{e},.{}.)",
                    if forward { "T" } else { "F" }
                ));
                format!("#{o}")
            })
            .collect();
        let lp = self.add(&format!("EDGE_LOOP('',({}))", oriented.join(",")));
        let bound = self.add(&format!("FACE_OUTER_BOUND('',#{lp},.T.)"));
        self.add(&format!(
            "ADVANCED_FACE('{name}',(#{bound}),#{surface},.T.)"
        ))
    }

    fn solid(&mut self, faces: &[u64]) -> u64 {
        let list: Vec<String> = faces.iter().map(|f| format!("#{f}")).collect();
        let shell = self.add(&format!("CLOSED_SHELL('',({}))", list.join(",")));
        self.add(&format!("MANIFOLD_SOLID_BREP('',#{shell})"))
    }

    fn text(&self) -> String {
        let mut s = String::from(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION((''),'2;1');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n",
        );
        for (k, l) in self.lines.iter().enumerate() {
            s += &format!("#{} = {l};\n", k + 1);
        }
        s + "ENDSEC;\nEND-ISO-10303-21;\n"
    }
}

const MILLIMETRE: &str = "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )";

/// Pavé a × b × c, une face par direction (`x0`, `x1`… si `named`, sinon sans nom).
fn cuboid(s: &mut Step, [a, b, c]: [f64; 3], named: bool) -> Vec<u64> {
    // sommet k = x + 2y + 4z
    let at = |k: usize| {
        [
            a * (k & 1) as f64,
            b * ((k >> 1) & 1) as f64,
            c * ((k >> 2) & 1) as f64,
        ]
    };
    let vertices: Vec<u64> = (0..8).map(|k| s.vertex(at(k))).collect();
    let mut edges = HashMap::new();
    for i in 0..8usize {
        for bit in [1, 2, 4] {
            if i & bit == 0 {
                let j = i | bit;
                edges.insert((i, j), s.line(vertices[i], vertices[j], at(i), at(j)));
            }
        }
    }
    // cycles antihoraires vus de dehors
    let cycles = [
        ("x0", [0, 4, 6, 2], [-1.0, 0.0, 0.0]),
        ("x1", [1, 3, 7, 5], [1.0, 0.0, 0.0]),
        ("y0", [0, 1, 5, 4], [0.0, -1.0, 0.0]),
        ("y1", [2, 6, 7, 3], [0.0, 1.0, 0.0]),
        ("z0", [0, 2, 3, 1], [0.0, 0.0, -1.0]),
        ("z1", [4, 5, 7, 6], [0.0, 0.0, 1.0]),
    ];
    cycles
        .iter()
        .map(|(name, cycle, normal)| {
            let loop_edges: Vec<(u64, bool)> = (0..4)
                .map(|k| {
                    let (p, q) = (cycle[k], cycle[(k + 1) % 4]);
                    (edges[&(p.min(q), p.max(q))], p < q)
                })
                .collect();
            let plane = s.plane(at(cycle[0]), *normal);
            s.face(if named { name } else { "" }, &loop_edges, plane)
        })
        .collect()
}

fn total_length(b: &Brep) -> f64 {
    b.edges.iter().map(|e| e.curve.length()).sum()
}

/// V − E + F = 2 pour une pièce sans trou ni couture.
fn euler(b: &Brep) -> i64 {
    b.vertices.len() as i64 - b.edges.len() as i64 + b.faces.len() as i64
}

#[test]
fn pave_lu_avec_ses_longueurs_et_ses_noms() {
    let mut s = Step::default();
    let faces = cuboid(&mut s, [2.0, 3.0, 4.0], true);
    s.solid(&faces);
    s.add(MILLIMETRE);
    let (b, unit) = read_brep(&s.text()).unwrap();
    assert_eq!(unit, "mm");
    assert_eq!((b.vertices.len(), b.edges.len(), b.faces.len()), (8, 12, 6));
    assert_eq!(euler(&b), 2);
    assert!((total_length(&b) - 4.0 * 9.0).abs() < 1e-12);
    let names: Vec<&str> = b.faces.iter().map(|f| f.name.as_str()).collect();
    assert_eq!(names, ["x0", "x1", "y0", "y1", "z0", "z1"]);
    assert!(b
        .faces
        .iter()
        .all(|f| f.loops.len() == 1 && f.loops[0].len() == 4));
    // chaque contour est fermé : la fin d'une arête est le début de la suivante
    for f in &b.faces {
        let lp = &f.loops[0];
        for k in 0..lp.len() {
            let ends = |(e, rev): (usize, bool)| {
                let edge = &b.edges[e];
                if rev {
                    (edge.end, edge.start)
                } else {
                    (edge.start, edge.end)
                }
            };
            assert_eq!(ends(lp[k]).1, ends(lp[(k + 1) % lp.len()]).0);
        }
    }
    // faces sans nom : face_1, face_2…
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], false);
    s.solid(&faces);
    let (b, unit) = read_brep(&s.text()).unwrap();
    assert_eq!(b.faces[0].name, "face_1");
    assert_eq!(b.faces[5].name, "face_6");
    assert_eq!(unit, "inconnue");
}

/// Quart de cylindre de rayon 2 et de hauteur 1 : l'arc du bas en B-spline rationnelle
/// (instance complexe, comme l'écrivent les modeleurs), celui du haut en CIRCLE.
fn quarter_cylinder(s: &mut Step, top_radius: f64) {
    let p = |x: f64, y: f64, z: f64| [x, y, z];
    let (o0, a0, b0) = (p(0., 0., 0.), p(2., 0., 0.), p(0., 2., 0.));
    let (o1, a1, b1) = (p(0., 0., 1.), p(2., 0., 1.), p(0., 2., 1.));
    let [vo0, va0, vb0, vo1, va1, vb1] = [o0, a0, b0, o1, a1, b1].map(|q| s.vertex(q));
    let oa0 = s.line(vo0, va0, o0, a0);
    let bo0 = s.line(vb0, vo0, b0, o0);
    let oa1 = s.line(vo1, va1, o1, a1);
    let bo1 = s.line(vb1, vo1, b1, o1);
    let oo = s.line(vo0, vo1, o0, o1);
    let aa = s.line(va0, va1, a0, a1);
    let bb = s.line(vb0, vb1, b0, b1);
    let control = [a0, p(2., 2., 0.), b0].map(|q| s.point(q));
    let spline = s.add(&format!(
        "( BOUNDED_CURVE() B_SPLINE_CURVE(2,(#{},#{},#{}),.CIRCULAR_ARC.,.F.,.F.) B_SPLINE_CURVE_WITH_KNOTS((3,3),(0.,1.),.PIECEWISE_BEZIER_KNOTS.) CURVE() GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_CURVE((1.,{FRAC_1_SQRT_2:?},1.)) REPRESENTATION_ITEM('') )",
        control[0], control[1], control[2]
    ));
    let ab0 = s.edge(va0, vb0, spline);
    let placement = s.placement(o1);
    let circle = s.add(&format!("CIRCLE('',#{placement},{top_radius:?})"));
    let ab1 = s.edge(va1, vb1, circle);
    let planes = [
        s.plane(o0, [0., 0., -1.]),
        s.plane(o1, [0., 0., 1.]),
        s.plane(o0, [0., -1., 0.]),
        s.plane(o0, [-1., 0., 0.]),
    ];
    let cylinder = s.cylinder(2.0);
    let faces = [
        s.face(
            "dessous",
            &[(bo0, false), (ab0, false), (oa0, false)],
            planes[0],
        ),
        s.face(
            "dessus",
            &[(oa1, true), (ab1, true), (bo1, true)],
            planes[1],
        ),
        s.face(
            "sym_y",
            &[(oa0, true), (aa, true), (oa1, false), (oo, false)],
            planes[2],
        ),
        s.face(
            "sym_x",
            &[(oo, true), (bo1, false), (bb, false), (bo0, true)],
            planes[3],
        ),
        s.face(
            "courbe",
            &[(ab0, true), (bb, true), (ab1, false), (aa, false)],
            cylinder,
        ),
    ];
    s.solid(&faces);
}

#[test]
fn quart_de_cylindre_arcs_exacts() {
    let mut s = Step::default();
    quarter_cylinder(&mut s, 2.0);
    let (b, _) = read_brep(&s.text()).unwrap();
    assert_eq!((b.vertices.len(), b.edges.len(), b.faces.len()), (6, 9, 5));
    assert_eq!(euler(&b), 2);
    // 4 segments de 2, 3 verticales de 1, deux quarts de cercle de rayon 2
    assert!(
        (total_length(&b) - (11.0 + 2.0 * PI)).abs() < 1e-9,
        "{}",
        total_length(&b)
    );
    let arcs: Vec<&EdgeCurve> = b
        .edges
        .iter()
        .map(|e| &e.curve)
        .filter(|c| !matches!(c, EdgeCurve::Line { .. }))
        .collect();
    assert_eq!(arcs.len(), 2);
    assert!(matches!(arcs[0], EdgeCurve::BSpline { .. }));
    assert!(matches!(arcs[1], EdgeCurve::Arc { .. }));
    for arc in arcs {
        assert!((arc.length() - PI).abs() < 1e-10);
        // l'arc va du sommet en x au sommet en y, et passe par 45°
        let [m0, m1, _] = arc.point(0.5);
        assert!((m0 - 2.0 * FRAC_1_SQRT_2).abs() < 1e-12 && (m1 - m0).abs() < 1e-12);
    }
}

/// Cylindre de rayon r et de hauteur h comme l'écrit OCCT : un sommet par cercle, des courbes
/// fermées et une couture (SEAM_CURVE) parcourue deux fois par la face latérale.
fn seamed_cylinder(s: &mut Step, r: f64, h: f64) {
    let (vb, vt) = (s.vertex([r, 0.0, 0.0]), s.vertex([r, 0.0, h]));
    let circle = |s: &mut Step, z: f64, v: u64| {
        let placement = s.placement([0.0, 0.0, z]);
        let c = s.add(&format!("CIRCLE('',#{placement},{r:?})"));
        s.edge(v, v, c)
    };
    let bottom = circle(s, 0.0, vb);
    let top = circle(s, h, vt);
    let p = s.point([r, 0.0, 0.0]);
    let d = s.direction([0.0, 0.0, 1.0]);
    let vector = s.add(&format!("VECTOR('',#{d},1.)"));
    let line = s.add(&format!("LINE('',#{p},#{vector})"));
    let seam_curve = s.add(&format!("SEAM_CURVE('',#{line},(),.PCURVE_S1.)"));
    let seam = s.edge(vb, vt, seam_curve);
    let side = s.cylinder(r);
    let (under, over) = (
        s.plane([0.0; 3], [0.0, 0.0, -1.0]),
        s.plane([0.0, 0.0, h], [0.0, 0.0, 1.0]),
    );
    let faces = [
        s.face(
            "cote",
            &[(bottom, true), (seam, true), (top, false), (seam, false)],
            side,
        ),
        s.face("fond", &[(bottom, false)], under),
        s.face("haut", &[(top, true)], over),
    ];
    s.solid(&faces);
}

#[test]
fn cylindre_a_couture() {
    let mut s = Step::default();
    let (r, h) = (2.0, 3.0);
    seamed_cylinder(&mut s, r, h);
    s.add("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )");
    let (brep, unit) = read_brep(&s.text()).unwrap();
    assert_eq!(unit, "m");
    assert_eq!(
        (brep.vertices.len(), brep.edges.len(), brep.faces.len()),
        (2, 3, 3)
    );
    assert_eq!(brep.edges[0].start, brep.edges[0].end);
    assert!((total_length(&brep) - (4.0 * PI * r + h)).abs() < 1e-12);
    assert_eq!(
        brep.faces[0].surface,
        FaceSurface::Cylinder {
            origin: [0.0; 3],
            axis: [0.0, 0.0, 1.0],
            radius: r,
        }
    );
    assert!(matches!(
        brep.faces[2].surface,
        FaceSurface::Plane {
            normal: [0.0, 0.0, 1.0],
            ..
        }
    ));
}

/// Volume des tétraèdres de la pièce maillée à la taille h.
fn meshed_volume(brep: &Brep, h: f64) -> f64 {
    let options = TetOptions {
        size: Some(h),
        ..TetOptions::default()
    };
    let t = tetrahedralize_brep(brep, &options).unwrap();
    let p = |i: usize| [t.nodes[i].x, t.nodes[i].y, t.nodes[i].z];
    t.tetrahedra
        .iter()
        .map(|e| {
            let [a, b, c, d] = e.map(p);
            let (u, v, w) = (
                [b[0] - a[0], b[1] - a[1], b[2] - a[2]],
                [c[0] - a[0], c[1] - a[1], c[2] - a[2]],
                [d[0] - a[0], d[1] - a[1], d[2] - a[2]],
            );
            (u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]))
                / 6.0
        })
        .sum()
}

#[test]
fn pieces_lues_puis_maillees() {
    // quart de cylindre de rayon 2 et de hauteur 1 : π ; cylindre à couture : 12π. Sommets sur
    // le cylindre : volume par défaut, de moins que la flèche relative des cordes (h / r)² / 8
    let sag = |h: f64| (h / 2.0f64).powi(2) / 8.0;
    let mut s = Step::default();
    quarter_cylinder(&mut s, 2.0);
    let (b, _) = read_brep(&s.text()).unwrap();
    let e = (PI - meshed_volume(&b, 0.25)) / PI;
    assert!(e > 0.0 && e < sag(0.25), "quart de cylindre : écart {e:e}");
    let mut s = Step::default();
    seamed_cylinder(&mut s, 2.0, 3.0);
    let (b, _) = read_brep(&s.text()).unwrap();
    let e = (12.0 * PI - meshed_volume(&b, 0.5)) / (12.0 * PI);
    assert!(e > 0.0 && e < sag(0.5), "cylindre : écart {e:e}");
}

#[test]
fn tronc_de_cone_en_degres() {
    // rayon 1 en bas, 2 en haut, hauteur 1 : demi-angle 45°, volume 7π/3
    let mut s = Step::default();
    let (vb, vt) = (s.vertex([1.0, 0.0, 0.0]), s.vertex([2.0, 0.0, 1.0]));
    let circle = |s: &mut Step, r: f64, z: f64, v: u64| {
        let placement = s.placement([0.0, 0.0, z]);
        let c = s.add(&format!("CIRCLE('',#{placement},{r:?})"));
        s.edge(v, v, c)
    };
    let bottom = circle(&mut s, 1.0, 0.0, vb);
    let top = circle(&mut s, 2.0, 1.0, vt);
    let seam = s.line(vb, vt, [1.0, 0.0, 0.0], [2.0, 0.0, 1.0]);
    let frame = s.placement([0.0; 3]);
    let cone = s.add(&format!("CONICAL_SURFACE('',#{frame},1.,45.)"));
    let (under, over) = (
        s.plane([0.0; 3], [0.0, 0.0, -1.0]),
        s.plane([0.0, 0.0, 1.0], [0.0, 0.0, 1.0]),
    );
    let faces = [
        s.face(
            "cote",
            &[(bottom, true), (seam, true), (top, false), (seam, false)],
            cone,
        ),
        s.face("fond", &[(bottom, false)], under),
        s.face("haut", &[(top, true)], over),
    ];
    s.solid(&faces);
    let radian = s.add("( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) )");
    let measure = s.add(&format!(
        "PLANE_ANGLE_MEASURE_WITH_UNIT(PLANE_ANGLE_MEASURE(0.017453292519943295),#{radian})"
    ));
    let dims = s.add("DIMENSIONAL_EXPONENTS(0.,0.,0.,0.,0.,0.,0.)");
    let degree = s.add(&format!(
        "( CONVERSION_BASED_UNIT('DEGREE',#{measure}) NAMED_UNIT(#{dims}) PLANE_ANGLE_UNIT() )"
    ));
    // le contexte dit quelle unité d'angle vaut pour la géométrie : le degré, pas le radian
    s.add(&format!(
        "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{degree})) REPRESENTATION_CONTEXT('','') )"
    ));
    let (b, _) = read_brep(&s.text()).unwrap();
    let FaceSurface::Cone { half_angle, .. } = b.faces[0].surface else {
        panic!("{:?}", b.faces[0].surface);
    };
    assert!((half_angle - PI / 4.0).abs() < 1e-15);
    let exact = 7.0 * PI / 3.0;
    let e = (exact - meshed_volume(&b, 0.25)) / exact;
    assert!(e > 0.0 && e < 1e-2, "écart {e:e}");
    // unité d'angle illisible : le cône est refusé
    s.lines[measure as usize - 1] = format!("PLANE_ANGLE_MEASURE_WITH_UNIT($,#{radian})");
    assert!(message(&s.text()).contains("unité d'angle"));
}

#[test]
fn pouces() {
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], true);
    s.solid(&faces);
    let measure = s.add("LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#999)");
    let dims = s.add("DIMENSIONAL_EXPONENTS(1.,0.,0.,0.,0.,0.,0.)");
    s.add(&format!(
        "( CONVERSION_BASED_UNIT('INCH',#{measure}) LENGTH_UNIT() NAMED_UNIT(#{dims}) )"
    ));
    assert_eq!(read_brep(&s.text()).unwrap().1, "inch");
    // une unité convertie au nom d'unité SI, comme l'écrit Onshape
    let last = s.lines.len() - 1;
    s.lines[last] = s.lines[last].replace("'INCH'", "'METRE'");
    assert_eq!(read_brep(&s.text()).unwrap().1, "m");
}

#[test]
fn sphere_et_tore() {
    // le sens de la face dit si la sphère est bombée : celle du cylindre arrondi l'est, l'intérieur
    // d'une coque ne l'est pas
    let mut s = Step::default();
    seamed_cylinder(&mut s, 2.0, 3.0);
    let top = s
        .lines
        .iter()
        .rposition(|l| l.starts_with("PLANE("))
        .unwrap();
    let frame = s.placement([0.0, 0.0, 3.0]);
    s.lines[top] = format!("SPHERICAL_SURFACE('',#{frame},2.)");
    let (b, _) = read_brep(&s.text()).unwrap();
    assert_eq!(
        b.faces[2].surface,
        FaceSurface::Sphere {
            center: [0.0, 0.0, 3.0],
            radius: 2.0,
            convex: true,
        }
    );
    // le bout arrondi se maille : 12π + 16π/3
    let exact = 12.0 * PI + 16.0 * PI / 3.0;
    let e = (exact - meshed_volume(&b, 0.5)) / exact;
    assert!(
        e > 0.0 && e < 2.0 * (0.5f64 / 2.0).powi(2) / 8.0,
        "écart {e:e}"
    );
    let face = s
        .lines
        .iter()
        .rposition(|l| l.starts_with("ADVANCED_FACE('haut'"))
        .unwrap();
    s.lines[face] = s.lines[face].replace(".T.)", ".F.)");
    let (b, _) = read_brep(&s.text()).unwrap();
    assert!(matches!(
        b.faces[2].surface,
        FaceSurface::Sphere { convex: false, .. }
    ));
    // tore : la lecture prend ses rayons (la géométrie est contrôlée au maillage)
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], true);
    s.solid(&faces);
    let plane = s
        .lines
        .iter()
        .position(|l| l.starts_with("PLANE("))
        .unwrap();
    s.lines[plane] = s.lines[plane]
        .replace("PLANE(", "TOROIDAL_SURFACE(")
        .replace(")", ",5.,1.5)");
    let (b, _) = read_brep(&s.text()).unwrap();
    assert!(matches!(
        b.faces[0].surface,
        FaceSurface::Torus {
            major: 5.0,
            minor: 1.5,
            ..
        }
    ));
}

fn message(text: &str) -> String {
    read_brep(text).unwrap_err().to_string()
}

#[test]
fn refus_avec_message() {
    // deux solides
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], true);
    let solid = s.solid(&faces);
    s.add(&format!("MANIFOLD_SOLID_BREP('',#{})", solid - 1));
    let m = message(&s.text());
    assert!(
        m.starts_with("STEP : 2 solides") && m.contains("une seule pièce"),
        "{m}"
    );
    // aucun solide
    let mut s = Step::default();
    cuboid(&mut s, [1.0; 3], true);
    assert!(message(&s.text()).contains("aucun solide"));
    // cavités
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], true);
    let solid = s.solid(&faces);
    s.lines[solid as usize - 1] = format!("BREP_WITH_VOIDS('',#{},())", solid - 1);
    assert!(message(&s.text()).contains("cavités"));
    // coque ouverte : une face en moins
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], true);
    s.solid(&faces[..5]);
    let m = message(&s.text());
    assert!(
        m.contains("parcourue 1 fois dans un sens et 0 dans l'autre")
            || m.contains("parcourue 0 fois dans un sens et 1"),
        "{m}"
    );
    // courbe non prise en charge
    let mut s = Step::default();
    let faces = cuboid(&mut s, [1.0; 3], true);
    s.solid(&faces);
    let line = s.lines.iter().position(|l| l.starts_with("LINE(")).unwrap();
    s.lines[line] = "OFFSET_CURVE_3D('',#1,1.,.F.,#2)".to_string();
    let m = message(&s.text());
    assert!(
        m.contains(&format!("ligne {} : courbe OFFSET_CURVE_3D", line + 7)),
        "{m}"
    );
    // surface non prise en charge
    let mut s = Step::default();
    quarter_cylinder(&mut s, 2.0);
    let cylinder = s
        .lines
        .iter()
        .position(|l| l.starts_with("CYLINDRICAL_SURFACE("))
        .unwrap();
    s.lines[cylinder] =
        s.lines[cylinder].replace("CYLINDRICAL_SURFACE", "SURFACE_OF_LINEAR_EXTRUSION");
    let m = message(&s.text());
    assert!(
        m.contains(&format!(
            "ligne {} : face sur une surface SURFACE_OF_LINEAR_EXTRUSION",
            cylinder + 7
        )),
        "{m}"
    );
    // sommet hors de son cercle
    let mut s = Step::default();
    quarter_cylinder(&mut s, 2.5);
    assert!(message(&s.text()).contains("n'est pas sur son cercle"));
    // référence absente : la coque cite une face qui n'existe pas
    let mut s = Step::default();
    s.solid(&[]);
    let shell = s.lines.len() - 2;
    s.lines[shell] = "CLOSED_SHELL('',(#9999))".to_string();
    let m = message(&s.text());
    assert!(m.contains("#9999 est citée"), "{m}");
}
