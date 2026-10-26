//! Prédicats géométriques en f64. Un résultat plus petit que l'erreur d'arrondi possible est
//! ramené à 0 : le mailleur traite alors le cas comme aligné ou cocirculaire au lieu de se fier à
//! un signe faux. Pas d'arithmétique exacte (Shewchuk) : un cas vraiment ambigu finit en erreur
//! `MeshingFailed`, jamais en triangle retourné.

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
}
