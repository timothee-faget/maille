// Mesh coloré : le GPU interpole la valeur du champ entre les sommets, puis chaque pixel lit sa
// couleur dans la LUT. Interpoler la valeur (et non la couleur) fait suivre aux bandes les isovaleurs.

struct Uniforms {
    // clip = scale * position + offset
    scale: vec2<f32>,
    offset: vec2<f32>,
    // (min, max) du champ affiché
    range: vec2<f32>,
}

@group(0) @binding(0) var<uniform> u: Uniforms;
@group(0) @binding(1) var lut: texture_1d<f32>;
@group(0) @binding(2) var lut_sampler: sampler;

struct VertexOut {
    @builtin(position) clip: vec4<f32>,
    @location(0) value: f32,
}

@vertex
fn vs_main(@location(0) position: vec2<f32>, @location(1) value: f32) -> VertexOut {
    var out: VertexOut;
    out.clip = vec4<f32>(position * u.scale + u.offset, 0.0, 1.0);
    out.value = value;
    return out;
}

// Légende : positions déjà en coordonnées de clip, ni zoom ni déplacement.
@vertex
fn vs_screen(@location(0) position: vec2<f32>, @location(1) value: f32) -> VertexOut {
    var out: VertexOut;
    out.clip = vec4<f32>(position, 0.0, 1.0);
    out.value = value;
    return out;
}

@fragment
fn fs_fill(in: VertexOut) -> @location(0) vec4<f32> {
    // champ constant : intervalle nul, on évite la division par zéro
    let span = max(u.range.y - u.range.x, 1e-30);
    let t = clamp((in.value - u.range.x) / span, 0.0, 1.0);
    return textureSample(lut, lut_sampler, t);
}

@fragment
fn fs_edge(in: VertexOut) -> @location(0) vec4<f32> {
    return vec4<f32>(0.15, 0.15, 0.15, 1.0);
}
