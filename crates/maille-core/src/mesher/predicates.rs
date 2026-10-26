//! Prédicats géométriques en f64. Un résultat plus petit que l'erreur d'arrondi possible est
//! ramené à 0 : le mailleur traite alors le cas comme aligné ou cocirculaire au lieu de se fier à
//! un signe faux. Pas d'arithmétique exacte (Shewchuk) : un cas vraiment ambigu finit en erreur
//! `MeshingFailed`, jamais en triangle ou tétraèdre retourné.

/// Marge relative sur l'erreur d'arrondi, large devant les bornes de Shewchuk (~10⁻¹⁵).
const EPS: f64 = 1e-12;

/// Deux fois l'aire signée du triangle abc : positive si a, b, c tournent dans le sens
/// antihoraire, négative dans le sens horaire, 0 si alignés (à l'arrondi près).
pub fn orient2d(a: [f64; 2], b: [f64; 2], c: [f64; 2]) -> f64 {
    let (acx, acy) = (a[0] - c[0], a[1] - c[1]);
    let (bcx, bcy) = (b[0] - c[0], b[1] - c[1]);
    let (l, r) = (acx * bcy, acy * bcx);
    filter(l - r, l.abs() + r.abs())
}

/// Positif si d est strictement dans le cercle circonscrit du triangle abc (antihoraire),
/// négatif s'il est dehors, 0 s'il est dessus (à l'arrondi près).
pub fn incircle(a: [f64; 2], b: [f64; 2], c: [f64; 2], d: [f64; 2]) -> f64 {
    let (adx, ady) = (a[0] - d[0], a[1] - d[1]);
    let (bdx, bdy) = (b[0] - d[0], b[1] - d[1]);
    let (cdx, cdy) = (c[0] - d[0], c[1] - d[1]);
    let (alift, blift, clift) = (
        adx * adx + ady * ady,
        bdx * bdx + bdy * bdy,
        cdx * cdx + cdy * cdy,
    );
    let (bc, cb) = (bdx * cdy, cdx * bdy);
    let (ca, ac) = (cdx * ady, adx * cdy);
    let (ab, ba) = (adx * bdy, bdx * ady);
    let det = alift * (bc - cb) + blift * (ca - ac) + clift * (ab - ba);
    let permanent = alift * (bc.abs() + cb.abs())
        + blift * (ca.abs() + ac.abs())
        + clift * (ab.abs() + ba.abs());
    filter(det, permanent)
}

/// Deux fois six fois le volume signé du tétraèdre abcd : positif si abcd est un tétraèdre
/// positif (comme `Tet4`, d est du côté de la normale (b − a) × (c − a)), négatif s'il est
/// inversé, 0 si les quatre points sont coplanaires (à l'arrondi près).
pub fn orient3d(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3]) -> f64 {
    let (det, permanent) = det3(sub(b, a), sub(c, a), sub(d, a));
    filter(det, permanent)
}

/// Positif si e est strictement dans la sphère circonscrite du tétraèdre positif abcd
/// (`orient3d(a, b, c, d) > 0`), négatif s'il est dehors, 0 s'il est dessus (à l'arrondi près).
pub fn insphere(a: [f64; 3], b: [f64; 3], c: [f64; 3], d: [f64; 3], e: [f64; 3]) -> f64 {
    // déterminant 4 × 4 des lignes (p − e, |p − e|²), développé selon la dernière colonne
    let rows = [a, b, c, d].map(|p| {
        let v = sub(p, e);
        (v, v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    });
    let [(va, la), (vb, lb), (vc, lc), (vd, ld)] = rows;
    let minors = [
        det3(vb, vc, vd),
        det3(va, vc, vd),
        det3(va, vb, vd),
        det3(va, vb, vc),
    ];
    let lifts = [la, lb, lc, ld];
    let (mut det, mut permanent) = (0.0, 0.0);
    for (i, ((m, m_perm), l)) in minors.into_iter().zip(lifts).enumerate() {
        let sign = if i % 2 == 0 { -1.0 } else { 1.0 };
        det += sign * l * m;
        permanent += l * m_perm;
    }
    // le déterminant est positif pour e dehors quand abcd est positif : on inverse
    filter(-det, permanent)
}

/// Distance entre deux points de l'espace.
pub fn distance3(p: [f64; 3], q: [f64; 3]) -> f64 {
    let v = sub(p, q);
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

fn sub(p: [f64; 3], q: [f64; 3]) -> [f64; 3] {
    [p[0] - q[0], p[1] - q[1], p[2] - q[2]]
}

/// Déterminant 3 × 3 des lignes u, v, w, et son « permanent » (somme des modules des termes).
fn det3(u: [f64; 3], v: [f64; 3], w: [f64; 3]) -> (f64, f64) {
    let terms = [
        u[0] * (v[1] * w[2] - v[2] * w[1]),
        -u[1] * (v[0] * w[2] - v[2] * w[0]),
        u[2] * (v[0] * w[1] - v[1] * w[0]),
    ];
    let abs = |a: f64, b: f64| a.abs() + b.abs();
    let perms = [
        u[0].abs() * abs(v[1] * w[2], v[2] * w[1]),
        u[1].abs() * abs(v[0] * w[2], v[2] * w[0]),
        u[2].abs() * abs(v[0] * w[1], v[1] * w[0]),
    ];
    (terms.iter().sum(), perms.iter().sum())
}

/// Distance entre deux points du plan.
pub fn distance(p: [f64; 2], q: [f64; 2]) -> f64 {
    (p[0] - q[0]).hypot(p[1] - q[1])
}

/// `value`, ou 0 si son module est sous l'erreur d'arrondi possible (`permanent` : somme des
/// modules des termes du calcul).
fn filter(value: f64, permanent: f64) -> f64 {
    if value.abs() <= EPS * permanent {
        0.0
    } else {
        value
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn orientation() {
        assert!(orient2d([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]) > 0.0);
        assert!(orient2d([0.0, 0.0], [0.0, 1.0], [1.0, 0.0]) < 0.0);
        assert!((orient2d([0.0, 0.0], [2.0, 0.0], [0.0, 3.0]) - 6.0).abs() < 1e-12);
    }

    #[test]
    fn aligned_points_give_zero() {
        assert_eq!(orient2d([0.0, 0.0], [1.0, 1.0], [2.0, 2.0]), 0.0);
        // alignés seulement à l'arrondi près : 0.1 et 0.3 ne sont pas représentables
        assert_eq!(orient2d([0.1, 0.1], [0.2, 0.2], [0.3, 0.3]), 0.0);
        // loin de l'origine, la marge suit l'ordre de grandeur des coordonnées
        assert_eq!(
            orient2d([1e6, 1e6], [1e6 + 0.1, 1e6 + 0.1], [1e6 + 0.3, 1e6 + 0.3]),
            0.0
        );
    }

    #[test]
    fn inside_outside_circle() {
        let (a, b, c) = ([0.0, 0.0], [1.0, 0.0], [0.0, 1.0]); // cercle de centre (½, ½)
        assert!(incircle(a, b, c, [0.5, 0.5]) > 0.0);
        assert!(incircle(a, b, c, [2.0, 2.0]) < 0.0);
        // le point (1, 1) est sur le cercle, y compris à l'arrondi près pour des coordonnées décimales
        assert_eq!(incircle(a, b, c, [1.0, 1.0]), 0.0);
        assert_eq!(
            incircle([0.0, 0.0], [0.1, 0.0], [0.0, 0.1], [0.1, 0.1]),
            0.0
        );
    }

    const O: [f64; 3] = [0.0, 0.0, 0.0];
    const X: [f64; 3] = [1.0, 0.0, 0.0];
    const Y: [f64; 3] = [0.0, 1.0, 0.0];
    const Z: [f64; 3] = [0.0, 0.0, 1.0];

    #[test]
    fn orientation_3d() {
        assert!((orient3d(O, X, Y, Z) - 1.0).abs() < 1e-12);
        assert!(orient3d(O, Y, X, Z) < 0.0);
        assert!(
            (orient3d(O, [2.0, 0.0, 0.0], [0.0, 3.0, 0.0], [0.0, 0.0, 4.0]) - 24.0).abs() < 1e-12
        );
        // permutation paire : même signe, impaire : signe opposé
        assert!(orient3d(X, O, Z, Y) > 0.0);
        assert!(orient3d(X, Y, Z, O) < 0.0);
    }

    #[test]
    fn coplanar_points_give_zero() {
        assert_eq!(orient3d(O, X, Y, [1.0, 1.0, 0.0]), 0.0);
        // coplanaires à l'arrondi près : 0.1 et 0.3 ne sont pas représentables
        assert_eq!(
            orient3d(
                [0.1, 0.1, 0.1],
                [0.2, 0.2, 0.2],
                [0.3, 0.3, 0.3],
                [0.1, 0.3, 0.2]
            ),
            0.0
        );
        assert_eq!(
            orient3d(
                [1e6, 1e6, 1e6],
                [1e6 + 0.1, 1e6, 1e6],
                [1e6, 1e6 + 0.1, 1e6],
                [1e6 + 0.3, 1e6 + 0.3, 1e6]
            ),
            0.0
        );
    }

    #[test]
    fn inside_outside_sphere() {
        // sphère de centre (½, ½, ½) et de rayon √3 / 2 : passe par les huit sommets du cube unité
        assert!(orient3d(O, X, Y, Z) > 0.0);
        assert!(insphere(O, X, Y, Z, [0.5, 0.5, 0.5]) > 0.0);
        assert!(insphere(O, X, Y, Z, [3.0, 3.0, 3.0]) < 0.0);
        assert!(insphere(O, X, Y, Z, [-1.0, 0.0, 0.0]) < 0.0);
        assert_eq!(insphere(O, X, Y, Z, [1.0, 1.0, 1.0]), 0.0);
        assert_eq!(insphere(O, X, Y, Z, [1.0, 1.0, 0.0]), 0.0);
        // un sommet est sur la sphère
        assert_eq!(insphere(O, X, Y, Z, Z), 0.0);
        // décimaux : (0.1, 0.1, 0.1) est cosphérique avec le tétraèdre de côté 0.1
        assert_eq!(
            insphere(
                O,
                [0.1, 0.0, 0.0],
                [0.0, 0.1, 0.0],
                [0.0, 0.0, 0.1],
                [0.1, 0.1, 0.1]
            ),
            0.0
        );
    }

    #[test]
    fn insphere_matches_distance_to_the_circumcentre() {
        let (a, b, c, d) = (
            [0.3, -0.2, 0.1],
            [2.1, 0.4, -0.3],
            [0.5, 1.7, 0.2],
            [0.6, 0.5, 1.9],
        );
        assert!(orient3d(a, b, c, d) > 0.0);
        // centre de la sphère par résolution de 2 (p − a) · o = |p|² − |a|²
        let rows = [b, c, d].map(|p| {
            let v = sub(p, a);
            (v, (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / 2.0)
        });
        let det = det3(rows[0].0, rows[1].0, rows[2].0).0;
        let col = |k: usize| {
            let mut m = [rows[0].0, rows[1].0, rows[2].0];
            for (row, r) in m.iter_mut().zip(&rows) {
                row[k] = r.1;
            }
            det3(m[0], m[1], m[2]).0 / det
        };
        let centre = [a[0] + col(0), a[1] + col(1), a[2] + col(2)];
        let radius = distance3(centre, a);
        for p in [b, c, d] {
            assert!((distance3(centre, p) - radius).abs() < 1e-12);
        }
        for dir in [
            [1.0, 0.0, 0.0],
            [0.0, -1.0, 0.0],
            [0.6, 0.0, 0.8],
            [-0.48, 0.64, -0.6],
        ] {
            let at = |r: f64| {
                [
                    centre[0] + r * dir[0],
                    centre[1] + r * dir[1],
                    centre[2] + r * dir[2],
                ]
            };
            assert!(insphere(a, b, c, d, at(0.9 * radius)) > 0.0);
            assert!(insphere(a, b, c, d, at(1.1 * radius)) < 0.0);
        }
    }
}
