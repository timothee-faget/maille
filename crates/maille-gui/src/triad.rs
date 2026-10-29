//! Trièdre d'orientation : où pointent x, y et z à l'écran pour la caméra courante.

/// Un axe du trièdre vu par la caméra.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Axis {
    pub name: &'static str,
    /// Extrémité de l'axe à l'écran, de longueur 1 pour un axe dans le plan de l'écran : x vers la
    /// droite, y vers le bas (comme les pixels).
    pub screen: [f64; 2],
    /// Composante vers l'observateur : positive, l'axe pointe vers nous ; négative, il s'éloigne.
    pub toward: f64,
}

/// Les trois axes pour la base de la caméra `(droite, haut, œil)`, du plus lointain au plus proche :
/// dessinés dans cet ordre, les plus proches passent par-dessus.
pub fn axes(right: [f64; 3], up: [f64; 3], eye: [f64; 3]) -> [Axis; 3] {
    let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    let mut axes = [("X", 0), ("Y", 1), ("Z", 2)].map(|(name, k)| {
        let mut unit = [0.0; 3];
        unit[k] = 1.0;
        Axis {
            name,
            screen: [dot(unit, right), -dot(unit, up)],
            toward: dot(unit, eye),
        }
    });
    axes.sort_by(|a, b| a.toward.total_cmp(&b.toward));
    axes
}

#[cfg(test)]
mod tests {
    use super::*;

    fn find<'a>(axes: &'a [Axis; 3], name: &str) -> &'a Axis {
        axes.iter().find(|a| a.name == name).unwrap()
    }

    #[test]
    fn the_front_view_shows_x_to_the_right_y_up_and_z_toward_us() {
        // vue de face : on regarde vers −z, x à droite, y en haut
        let a = axes([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]);
        assert_eq!(find(&a, "X").screen, [1.0, 0.0]);
        assert_eq!(find(&a, "Y").screen, [0.0, -1.0]);
        let z = find(&a, "Z");
        assert_eq!(z.screen, [0.0, 0.0]);
        assert_eq!(z.toward, 1.0);
        assert_eq!(
            a[2].name, "Z",
            "l'axe le plus proche est dessiné en dernier"
        );
    }

    #[test]
    fn a_view_from_plus_x_points_x_at_us_and_z_to_the_left() {
        // œil en +x : droite = −z, haut = y
        let a = axes([0.0, 0.0, -1.0], [0.0, 1.0, 0.0], [1.0, 0.0, 0.0]);
        assert_eq!(find(&a, "X").toward, 1.0);
        assert_eq!(find(&a, "Z").screen, [-1.0, 0.0]);
    }

    #[test]
    fn an_axis_that_moves_away_comes_first() {
        let a = axes([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]);
        assert_eq!(a[0].name, "Z");
        assert!(a[0].toward < 0.0);
    }
}
