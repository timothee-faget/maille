//! Import STEP : topologie B-rep et courbes des arêtes. Pièces écrites à la main (pavé, quart de
//! cylindre à arcs en CIRCLE et en B-spline rationnelle, cylindre elliptique à couture comme chez
//! OCCT), validées par les longueurs exactes des arêtes, la formule d'Euler et l'usage des arêtes ;
//! unités et refus avec message.

use std::collections::HashMap;
use std::f64::consts::{FRAC_1_SQRT_2, PI};

use maille_io::step::{read_brep, Brep, EdgeCurve};

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
        let (o, z, x) = (
            self.point(origin),
            self.direction([0.0, 0.0, 1.0]),
            self.direction([1.0, 0.0, 0.0]),
        );
        self.add(&format!("AXIS2_PLACEMENT_3D('',#{o},#{z},#{x})"))
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

    /// Face nommée `name` : un contour extérieur d'arêtes (`true` : dans leur sens).
    fn face(&mut self, name: &str, edges: &[(u64, bool)]) -> u64 {
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
        let placement = self.placement([0.0; 3]);
        let plane = self.add(&format!("PLANE('',#{placement})"));
        self.add(&format!("ADVANCED_FACE('{name}',(#{bound}),#{plane},.T.)"))
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
        ("x0", [0, 4, 6, 2]),
        ("x1", [1, 3, 7, 5]),
        ("y0", [0, 1, 5, 4]),
        ("y1", [2, 6, 7, 3]),
        ("z0", [0, 2, 3, 1]),
        ("z1", [4, 5, 7, 6]),
    ];
    cycles
        .iter()
        .map(|(name, cycle)| {
            let loop_edges: Vec<(u64, bool)> = (0..4)
                .map(|k| {
                    let (p, q) = (cycle[k], cycle[(k + 1) % 4]);
                    (edges[&(p.min(q), p.max(q))], p < q)
                })
                .collect();
            s.face(if named { name } else { "" }, &loop_edges)
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
    let faces = [
        s.face("dessous", &[(bo0, false), (ab0, false), (oa0, false)]),
        s.face("dessus", &[(oa1, true), (ab1, true), (bo1, true)]),
        s.face(
            "sym_y",
            &[(oa0, true), (aa, true), (oa1, false), (oo, false)],
        ),
        s.face(
            "sym_x",
            &[(oo, true), (bo1, false), (bb, false), (bo0, true)],
        ),
        s.face(
            "courbe",
            &[(ab0, true), (bb, true), (ab1, false), (aa, false)],
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

#[test]
fn cylindre_elliptique_a_couture() {
    // comme OCCT : un sommet par cercle, des courbes fermées et une couture (SEAM_CURVE) parcourue
    // deux fois par la face latérale
    let mut s = Step::default();
    let (a, b, h) = (2.0, 1.0, 3.0);
    let (vb, vt) = (s.vertex([a, 0.0, 0.0]), s.vertex([a, 0.0, h]));
    let ellipse = |s: &mut Step, z: f64, v: u64| {
        let placement = s.placement([0.0, 0.0, z]);
        let e = s.add(&format!("ELLIPSE('',#{placement},{a:?},{b:?})"));
        s.edge(v, v, e)
    };
    let bottom = ellipse(&mut s, 0.0, vb);
    let top = ellipse(&mut s, h, vt);
    let p = s.point([a, 0.0, 0.0]);
    let d = s.direction([0.0, 0.0, 1.0]);
    let vector = s.add(&format!("VECTOR('',#{d},1.)"));
    let line = s.add(&format!("LINE('',#{p},#{vector})"));
    let seam_curve = s.add(&format!("SEAM_CURVE('',#{line},(),.PCURVE_S1.)"));
    let seam = s.edge(vb, vt, seam_curve);
    let faces = [
        s.face(
            "cote",
            &[(bottom, true), (seam, true), (top, false), (seam, false)],
        ),
        s.face("fond", &[(bottom, false)]),
        s.face("haut", &[(top, true)]),
    ];
    s.solid(&faces);
    s.add("( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )");
    let (brep, unit) = read_brep(&s.text()).unwrap();
    assert_eq!(unit, "m");
    assert_eq!(
        (brep.vertices.len(), brep.edges.len(), brep.faces.len()),
        (2, 3, 3)
    );
    assert_eq!(brep.edges[0].start, brep.edges[0].end);
    // périmètre de l'ellipse de demi-axes 2 et 1 : 9,688448220547675…
    let want = 2.0 * 9.688_448_220_547_675 + h;
    assert!(
        (total_length(&brep) - want).abs() < 1e-8,
        "{}",
        total_length(&brep)
    );
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
