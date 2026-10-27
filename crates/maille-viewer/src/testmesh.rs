//! Maillages de tétraèdres pour les tests : un cube de côté 2 centré à l'origine, maillé par le
//! mailleur 3D du cœur (Tet4, ou Tet10 par `to_quadratic`).

use maille_core::mesher::{tetrahedralize, Patch, Surface};

/// (positions, cellules) du cube [-1, 1]³ en Tet4.
pub fn cube() -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
    let t = surface_cube();
    let mesh = tetrahedralize(&t).unwrap();
    (
        mesh.nodes.iter().map(|n| [n.x, n.y, n.z]).collect(),
        mesh.tetrahedra.iter().map(|c| c.to_vec()).collect(),
    )
}

/// Même cube en Tet10.
pub fn cube_quadratic() -> (Vec<[f64; 3]>, Vec<Vec<usize>>) {
    let mesh = tetrahedralize(&surface_cube()).unwrap().to_quadratic();
    (
        mesh.nodes.iter().map(|n| [n.x, n.y, n.z]).collect(),
        mesh.tetrahedra.iter().map(|c| c.to_vec()).collect(),
    )
}

fn surface_cube() -> Surface {
    let vertices: Vec<[f64; 3]> = (0..8)
        .map(|k| {
            [
                if k & 1 == 0 { -1.0 } else { 1.0 },
                if k & 2 == 0 { -1.0 } else { 1.0 },
                if k & 4 == 0 { -1.0 } else { 1.0 },
            ]
        })
        .collect();
    // faces -x, +x, -y, +y, -z, +z, deux triangles chacune, normales sortantes
    let quads: [[usize; 4]; 6] = [
        [0, 4, 6, 2],
        [1, 3, 7, 5],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 2, 3, 1],
        [4, 5, 7, 6],
    ];
    let triangles: Vec<[usize; 3]> = quads
        .iter()
        .flat_map(|q| [[q[0], q[1], q[2]], [q[0], q[2], q[3]]])
        .collect();
    Surface {
        vertices,
        patches: vec![Patch {
            name: "peau".to_string(),
            triangles: (0..triangles.len()).collect(),
        }],
        triangles,
    }
}
