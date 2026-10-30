// Mesh coloré : le GPU interpole la valeur du champ entre les sommets, puis chaque pixel lit sa
// couleur dans la LUT. Interpoler la valeur (et non la couleur) fait suivre aux bandes les isovaleurs.

struct Uniforms {
    // clip = view_proj * position (caméra orthographique, profondeur dans [0, 1])
    view_proj: mat4x4<f32>,
    // (min, max) du champ affiché
    range: vec2<f32>,
    // direction de l'œil : de la cible vers l'observateur (xyz)
    eye: vec4<f32>,
    // couleur des arêtes
    edge: vec4<f32>,
}

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var lut: texture_1d<f32>;
@group(0) @binding(2) var lut_sampler: sampler;

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) value: f32,
    @location(1) world: vec3<f32>,
}

@vertex
fn vs_main(@location(0) position: vec3<f32>, @location(1) value: f32) -> VertexOut {
    var out: VertexOut;
    out.clip = u.view_proj * vec4<f32>(position, 1.0);
    out.value = value;
    out.world = position;
    return out;
}

// Segments du contour : position seule.
@vertex
fn vs_line(@location(0) position: vec3<f32>) -> VertexOut {
    var out: VertexOut;
    out.clip = u.view_proj * vec4<f32>(position, 1.0);
    out.value = 0.0;
    out.world = position;
    return out;
}

// Légende : positions déjà en coordonnées de clip, ni zoom ni déplacement.
@vertex
fn vs_screen(@location(0) position: vec2<f32>, @location(1) value: f32) -> VertexOut {
    var out: VertexOut;
    out.clip = vec4<f32>(position, 0.0, 1.0);
    out.value = value;
    out.world = vec3<f32>(0.0);
    return out;
}

fn field_color(value: f32) -> vec4<f32> {
    // champ constant : intervalle nul, on évite la division par zéro
    let span = max(u.range.y - u.range.x, 1e-30);
    let t = clamp((value - u.range.x) / span, 0.0, 1.0);
    return textureSample(lut, lut_sampler, t);
}

@fragment
fn fs_fill(in: VertexOut) -> @location(0) vec4<f32> {
    return field_color(in.value);
}

// Mesh : même couleur, assombrie quand la face est vue de biais. La normale vient des dérivées de
// la position, sans buffer de normales. Une face vue de face (tout le 2D) garde sa couleur exacte.
@fragment
fn fs_shaded(in: VertexOut) -> @location(0) vec4<f32> {
    let color = field_color(in.value);
    let n = cross(dpdx(in.world), dpdy(in.world));
    let len = length(n);
    var shade = 1.0;
    if (len > 0.0) {
        shade = 0.55 + 0.45 * abs(dot(n / len, u.eye.xyz));
    }
    return vec4<f32>(color.rgb * shade, color.a);
}

// Faces surlignées : une couleur par sommet, mélangée au Mesh selon son alpha.
struct HighlightOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) color: vec4<f32>,
}

@vertex
fn vs_highlight(@location(0) position: vec3<f32>, @location(1) color: vec4<f32>) -> HighlightOut {
    var out: HighlightOut;
    out.clip = u.view_proj * vec4<f32>(position, 1.0);
    out.color = color;
    return out;
}

@fragment
fn fs_highlight(in: HighlightOut) -> @location(0) vec4<f32> {
    return in.color;
}

@fragment
fn fs_edge(in: VertexOut) -> @location(0) vec4<f32> {
    return u.edge;
}
