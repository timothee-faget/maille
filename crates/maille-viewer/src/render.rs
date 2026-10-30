//! Contexte GPU (wgpu) : pipelines, buffers du Mesh, dessin. La cible est une fenêtre (surface) ou
//! une texture hors écran, dont on relit les pixels pour une capture.

use std::error::Error;
use std::sync::Arc;

use wgpu::util::DeviceExt;
use winit::event_loop::ActiveEventLoop;
use winit::window::Window;

use crate::camera::Camera;
use crate::colormap::Colormap;
use crate::legend::{Layout, FONT_SIZE, TITLE_FONT_SIZE};

/// Apparence du viewer : fond, texte de la légende et arêtes. Le fond est un gris (jamais du blanc pur),
/// choisi pour que les couleurs de la colormap restent lisibles.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Theme {
    #[default]
    Light,
    Dark,
}

impl Theme {
    /// Fond, en sRGB 8 bits.
    fn background_srgb(self) -> [u8; 3] {
        match self {
            Theme::Light => [204, 204, 204],
            Theme::Dark => [37, 39, 43],
        }
    }

    /// Fond pour `LoadOp::Clear` : wgpu attend du linéaire, même sur une cible sRGB.
    fn background(self) -> wgpu::Color {
        let [r, g, b] = self.background_srgb().map(|c| {
            let c = f64::from(c) / 255.0;
            if c <= 0.04045 {
                c / 12.92
            } else {
                ((c + 0.055) / 1.055).powf(2.4)
            }
        });
        wgpu::Color { r, g, b, a: 1.0 }
    }

    fn text(self) -> glyphon::Color {
        match self {
            Theme::Light => glyphon::Color::rgb(38, 38, 38),
            Theme::Dark => glyphon::Color::rgb(225, 225, 225),
        }
    }

    /// Couleur des arêtes et du cadre de la légende (RGBA, sRGB).
    fn edge(self) -> [f32; 4] {
        match self {
            Theme::Light => [0.15, 0.15, 0.15, 1.0],
            Theme::Dark => [0.85, 0.85, 0.85, 1.0],
        }
    }
}

const DEPTH_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;
/// Format des captures hors écran : sRGB, comme la fenêtre, pour des couleurs identiques.
const OFFSCREEN_FORMAT: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;

enum Target {
    Window {
        window: Arc<Window>,
        surface: wgpu::Surface<'static>,
        config: wgpu::SurfaceConfiguration,
    },
    /// Texture hors écran : capture PNG, ou image montrée par un autre client (l'interface graphique)
    /// qui partage le device. `scale` : pixels physiques par pixel logique de cet écran.
    Offscreen {
        texture: wgpu::Texture,
        width: u32,
        height: u32,
        scale: f64,
    },
}

pub struct Renderer {
    target: Target,
    device: wgpu::Device,
    queue: wgpu::Queue,
    fill_pipeline: wgpu::RenderPipeline,
    edge_pipeline: wgpu::RenderPipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    sampler: wgpu::Sampler,
    /// Uniforms, LUT de la colormap et son sampler.
    bind_group: wgpu::BindGroup,
    /// Matrice de la caméra, (min, max) du champ et direction de l'œil : voir `Uniforms` dans
    /// shader.wgsl.
    uniforms: wgpu::Buffer,
    depth: wgpu::TextureView,
    mesh: Option<MeshBuffers>,
    camera: [[f32; 4]; 4],
    eye: [f32; 3],
    theme: Theme,
    range: (f64, f64),
    /// Barre de couleurs (`vs_screen` + `fs_fill`) et ses traits (`vs_screen` + `fs_edge`).
    legend_fill_pipeline: wgpu::RenderPipeline,
    legend_line_pipeline: wgpu::RenderPipeline,
    legend: Option<LegendBuffers>,
    /// Triangles surlignés par-dessus le Mesh : positions, couleurs (RGBA) et nombre de sommets.
    highlight_pipeline: wgpu::RenderPipeline,
    highlight: Option<(wgpu::Buffer, wgpu::Buffer, u32)>,
    text: Text,
}

/// Ce que le Mesh envoie au GPU : sommets (position, valeur du champ), triangles indexés, et les
/// extrémités des segments du contour.
struct MeshBuffers {
    positions: wgpu::Buffer,
    values: wgpu::Buffer,
    triangles: wgpu::Buffer,
    triangle_indices: u32,
    lines: wgpu::Buffer,
    line_vertices: u32,
}

/// Géométrie de la légende, en coordonnées de clip : les 4 coins de la barre (valeurs min en bas,
/// max en haut, colorées par le même shader que le Mesh), puis 2 sommets par graduation.
struct LegendBuffers {
    positions: wgpu::Buffer,
    values: wgpu::Buffer,
    bar: wgpu::Buffer,
    lines: wgpu::Buffer,
    line_indices: u32,
}

/// Texte de la légende (glyphon) : polices, cache de glyphes, et un bloc de texte par position.
struct Text {
    fonts: glyphon::FontSystem,
    swash_cache: glyphon::SwashCache,
    viewport: glyphon::Viewport,
    atlas: glyphon::TextAtlas,
    renderer: glyphon::TextRenderer,
    blocks: Vec<(glyphon::Buffer, [f32; 2])>,
}

fn depth_state(write: bool, compare: wgpu::CompareFunction, bias: bool) -> wgpu::DepthStencilState {
    wgpu::DepthStencilState {
        format: DEPTH_FORMAT,
        depth_write_enabled: Some(write),
        depth_compare: Some(compare),
        stencil: Default::default(),
        // le remplissage recule un peu, pour que le contour tracé sur la même surface reste visible
        bias: if bias {
            wgpu::DepthBiasState {
                constant: 2,
                slope_scale: 1.0,
                clamp: 0.0,
            }
        } else {
            Default::default()
        },
    }
}

fn depth_view(device: &wgpu::Device, width: u32, height: u32) -> wgpu::TextureView {
    device
        .create_texture(&wgpu::TextureDescriptor {
            label: None,
            size: wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D2,
            format: DEPTH_FORMAT,
            usage: wgpu::TextureUsages::RENDER_ATTACHMENT,
            view_formats: &[],
        })
        .create_view(&Default::default())
}

fn offscreen_texture(device: &wgpu::Device, width: u32, height: u32) -> wgpu::Texture {
    device.create_texture(&wgpu::TextureDescriptor {
        label: None,
        size: wgpu::Extent3d {
            width,
            height,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format: OFFSCREEN_FORMAT,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT
            | wgpu::TextureUsages::COPY_SRC
            | wgpu::TextureUsages::TEXTURE_BINDING,
        view_formats: &[],
    })
}

impl Renderer {
    /// Prépare le GPU pour une fenêtre ; la géométrie, le champ et la caméra sont ensuite fixés
    /// par `set_geometry`, `set_range` et `set_camera`.
    pub fn new_window(
        event_loop: &ActiveEventLoop,
        window: Arc<Window>,
        colormap: Colormap,
    ) -> Result<Self, Box<dyn Error>> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_with_display_handle(
            Box::new(event_loop.owned_display_handle()),
        ));
        let surface = instance.create_surface(window.clone())?;
        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            compatible_surface: Some(&surface),
            ..Default::default()
        }))
        .map_err(|e| format!("aucune carte graphique utilisable pour le viewer ({e})"))?;
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))?;
        let size = window.inner_size();
        let config = surface
            .get_default_config(&adapter, size.width.max(1), size.height.max(1))
            .ok_or("la carte graphique ne sait pas afficher dans cette fenêtre")?;
        surface.configure(&device, &config);
        let format = config.format;
        let target = Target::Window {
            window,
            surface,
            config,
        };
        Ok(Self::build(device, queue, format, target, colormap))
    }

    /// Prépare le GPU pour une image de `width` × `height` pixels, sans fenêtre.
    pub fn new_offscreen(
        width: u32,
        height: u32,
        colormap: Colormap,
    ) -> Result<Self, Box<dyn Error>> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor::new_without_display_handle());
        let adapter = pollster::block_on(instance.request_adapter(&Default::default()))
            .map_err(|e| format!("aucune carte graphique utilisable pour la capture ({e})"))?;
        let (device, queue) =
            pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor::default()))?;
        let (width, height) = (width.max(1), height.max(1));
        Ok(Self::new_shared(
            device, queue, width, height, 1.0, colormap,
        ))
    }

    /// Prépare le GPU pour une texture de `width` × `height` pixels sur un device qui appartient à un
    /// autre client (l'interface graphique) : `draw` la remplit, `texture` la donne. `scale` :
    /// pixels physiques par pixel logique, pour la taille du texte.
    pub fn new_shared(
        device: wgpu::Device,
        queue: wgpu::Queue,
        width: u32,
        height: u32,
        scale: f64,
        colormap: Colormap,
    ) -> Self {
        let (width, height) = (width.max(1), height.max(1));
        let target = Target::Offscreen {
            texture: offscreen_texture(&device, width, height),
            width,
            height,
            scale,
        };
        Self::build(device, queue, OFFSCREEN_FORMAT, target, colormap)
    }

    fn build(
        device: wgpu::Device,
        queue: wgpu::Queue,
        format: wgpu::TextureFormat,
        target: Target,
        colormap: Colormap,
    ) -> Self {
        let shader = device.create_shader_module(wgpu::include_wgsl!("shader.wgsl"));
        let bind_group_layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: None,
            entries: &[
                wgpu::BindGroupLayoutEntry {
                    binding: 0,
                    visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                    ty: wgpu::BindingType::Buffer {
                        ty: wgpu::BufferBindingType::Uniform,
                        has_dynamic_offset: false,
                        min_binding_size: None,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 1,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Texture {
                        sample_type: wgpu::TextureSampleType::Float { filterable: true },
                        view_dimension: wgpu::TextureViewDimension::D1,
                        multisampled: false,
                    },
                    count: None,
                },
                wgpu::BindGroupLayoutEntry {
                    binding: 2,
                    visibility: wgpu::ShaderStages::FRAGMENT,
                    ty: wgpu::BindingType::Sampler(wgpu::SamplerBindingType::Filtering),
                    count: None,
                },
            ],
        });
        let layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: None,
            bind_group_layouts: &[Some(&bind_group_layout)],
            immediate_size: 0,
        });
        let position3 = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<[f32; 3]>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![0 => Float32x3],
        };
        let position2 = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<[f32; 2]>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![0 => Float32x2],
        };
        let value = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<f32>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![1 => Float32],
        };
        let color = wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<[u8; 4]>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &wgpu::vertex_attr_array![1 => Unorm8x4],
        };
        let pipeline = |vertex,
                        buffers: &[Option<wgpu::VertexBufferLayout>],
                        topology,
                        fragment,
                        depth,
                        blend: Option<wgpu::BlendState>| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: None,
                layout: Some(&layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some(vertex),
                    compilation_options: Default::default(),
                    buffers,
                },
                // les deux faces se voient : la coupe montre l'intérieur de la peau
                primitive: wgpu::PrimitiveState {
                    topology,
                    cull_mode: None,
                    ..Default::default()
                },
                depth_stencil: Some(depth),
                multisample: Default::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fragment),
                    compilation_options: Default::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        blend,
                        ..format.into()
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let (triangle_list, line_list) = (
            wgpu::PrimitiveTopology::TriangleList,
            wgpu::PrimitiveTopology::LineList,
        );
        let mesh_buffers = [Some(position3.clone()), Some(value.clone())];
        let line_buffers = [Some(position3.clone())];
        let screen_buffers = [Some(position2), Some(value)];
        let lequal = wgpu::CompareFunction::LessEqual;
        let always = wgpu::CompareFunction::Always;
        let fill_pipeline = pipeline(
            "vs_main",
            &mesh_buffers,
            triangle_list,
            "fs_shaded",
            depth_state(true, lequal, true),
            None,
        );
        let edge_pipeline = pipeline(
            "vs_line",
            &line_buffers,
            line_list,
            "fs_edge",
            depth_state(true, lequal, false),
            None,
        );
        let highlight_pipeline = pipeline(
            "vs_highlight",
            &[Some(position3), Some(color)],
            triangle_list,
            "fs_highlight",
            depth_state(false, lequal, false),
            Some(wgpu::BlendState::ALPHA_BLENDING),
        );
        let legend_fill_pipeline = pipeline(
            "vs_screen",
            &screen_buffers,
            triangle_list,
            "fs_fill",
            depth_state(false, always, false),
            None,
        );
        let legend_line_pipeline = pipeline(
            "vs_screen",
            &screen_buffers,
            line_list,
            "fs_edge",
            depth_state(false, always, false),
            None,
        );

        let fonts = glyphon::FontSystem::new();
        let cache = glyphon::Cache::new(&device);
        let mut atlas = glyphon::TextAtlas::new(&device, &queue, &cache, format);
        let text = Text {
            // le texte passe dans la même passe que le Mesh : même format de profondeur, jamais testée
            renderer: glyphon::TextRenderer::new(
                &mut atlas,
                &device,
                Default::default(),
                Some(depth_state(false, always, false)),
            ),
            viewport: glyphon::Viewport::new(&device, &cache),
            swash_cache: glyphon::SwashCache::new(),
            blocks: Vec::new(),
            fonts,
            atlas,
        };

        let sampler = device.create_sampler(&wgpu::SamplerDescriptor::default());
        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            // matrice 4 × 4, (min, max) complété à 16 octets, direction de l'œil (vec4), couleur des arêtes (vec4)
            size: std::mem::size_of::<[f32; 28]>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let bind_group = lut_bind_group(
            &device,
            &queue,
            &bind_group_layout,
            &uniforms,
            &sampler,
            colormap,
        );
        let (width, height) = match &target {
            Target::Window { config, .. } => (config.width, config.height),
            Target::Offscreen { width, height, .. } => (*width, *height),
        };
        Self {
            depth: depth_view(&device, width, height),
            target,
            device,
            queue,
            fill_pipeline,
            edge_pipeline,
            bind_group_layout,
            sampler,
            bind_group,
            uniforms,
            mesh: None,
            camera: Camera::fit(([0.0; 3], [1.0; 3]), [1.0, 1.0], 0.0).matrix([1.0, 1.0]),
            eye: [0.0, 0.0, 1.0],
            theme: Theme::default(),
            range: (0.0, 1.0),
            legend_fill_pipeline,
            legend_line_pipeline,
            legend: None,
            highlight_pipeline,
            highlight: None,
            text,
        }
    }

    fn request_redraw(&self) {
        if let Target::Window { window, .. } = &self.target {
            window.request_redraw();
        }
    }

    /// Facteur d'échelle de l'écran : 2 sur un écran Retina (pixels physiques par pixel logique).
    pub fn scale_factor(&self) -> f64 {
        match &self.target {
            Target::Window { window, .. } => window.scale_factor(),
            Target::Offscreen { scale, .. } => *scale,
        }
    }

    /// Légende : barre de couleurs et graduations de `layout`, et `title` au-dessus (une ligne par
    /// `\n`). La barre prend le (min, max) du dernier `set_values`.
    pub fn set_legend(&mut self, layout: &Layout, title: &str) {
        let [w, h] = self.size();
        let clip = |x: f64, y: f64| [(2.0 * x / w - 1.0) as f32, (1.0 - 2.0 * y / h) as f32];
        let [left, top, right, bottom] = layout.bar;
        let (lo, hi) = (self.range.0 as f32, self.range.1 as f32);
        let mut positions = vec![
            clip(left, bottom),
            clip(right, bottom),
            clip(right, top),
            clip(left, top),
        ];
        let mut values = vec![lo, lo, hi, hi];
        let mut lines: Vec<u32> = vec![0, 1, 1, 2, 2, 3, 3, 0];
        for tick in &layout.ticks {
            lines.extend([positions.len() as u32, positions.len() as u32 + 1]);
            positions.extend([
                clip(right, tick.y),
                clip(right + layout.tick_length, tick.y),
            ]);
            values.extend([0.0, 0.0]);
        }
        let buffer = |contents: &[u8], usage| {
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: None,
                    contents,
                    usage,
                })
        };
        self.legend = Some(LegendBuffers {
            positions: buffer(bytemuck::cast_slice(&positions), wgpu::BufferUsages::VERTEX),
            values: buffer(bytemuck::cast_slice(&values), wgpu::BufferUsages::VERTEX),
            bar: buffer(
                bytemuck::cast_slice(&[0_u32, 1, 2, 0, 2, 3]),
                wgpu::BufferUsages::INDEX,
            ),
            lines: buffer(bytemuck::cast_slice(&lines), wgpu::BufferUsages::INDEX),
            line_indices: lines.len() as u32,
        });

        let scale = self.scale_factor();
        let mut blocks = vec![self.text_block(title, TITLE_FONT_SIZE * scale, layout.title)];
        for tick in &layout.ticks {
            let size = FONT_SIZE * scale;
            // étiquette centrée verticalement sur sa graduation
            let block = self.text_block(&tick.label, size, [layout.label_x, tick.y - 0.6 * size]);
            blocks.push(block);
        }
        self.text.blocks = blocks;
        self.request_redraw();
    }

    /// Bloc de texte de `size` pixels dont le coin haut gauche est en `at`.
    fn text_block(&mut self, text: &str, size: f64, at: [f64; 2]) -> (glyphon::Buffer, [f32; 2]) {
        let fonts = &mut self.text.fonts;
        let mut buffer =
            glyphon::Buffer::new(fonts, glyphon::Metrics::new(size as f32, 1.2 * size as f32));
        buffer.set_size(None, None);
        buffer.set_text(
            text,
            &glyphon::Attrs::new().family(glyphon::Family::SansSerif),
            glyphon::Shaping::Advanced,
            None,
        );
        buffer.shape_until_scroll(fonts, false);
        (buffer, [at[0] as f32, at[1] as f32])
    }

    pub fn set_colormap(&mut self, colormap: Colormap) {
        self.bind_group = lut_bind_group(
            &self.device,
            &self.queue,
            &self.bind_group_layout,
            &self.uniforms,
            &self.sampler,
            colormap,
        );
        self.request_redraw();
    }

    /// Taille de la zone de dessin, en pixels.
    pub fn size(&self) -> [f64; 2] {
        match &self.target {
            Target::Window { config, .. } => [f64::from(config.width), f64::from(config.height)],
            Target::Offscreen { width, height, .. } => [f64::from(*width), f64::from(*height)],
        }
    }

    /// Géométrie du Mesh : sommets (position et valeur du champ), triangles (indices) et segments
    /// du contour (deux extrémités chacun).
    pub fn set_geometry(
        &mut self,
        positions: &[[f32; 3]],
        values: &[f32],
        triangles: &[u32],
        lines: &[[f32; 3]],
    ) {
        let buffer = |contents: &[u8], usage| {
            self.device
                .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                    label: None,
                    contents,
                    usage,
                })
        };
        self.mesh = Some(MeshBuffers {
            positions: buffer(bytemuck::cast_slice(positions), wgpu::BufferUsages::VERTEX),
            values: buffer(bytemuck::cast_slice(values), wgpu::BufferUsages::VERTEX),
            triangles: buffer(bytemuck::cast_slice(triangles), wgpu::BufferUsages::INDEX),
            triangle_indices: triangles.len() as u32,
            lines: buffer(bytemuck::cast_slice(lines), wgpu::BufferUsages::VERTEX),
            line_vertices: lines.len() as u32,
        });
        self.request_redraw();
    }

    /// Triangles à surligner (3 sommets chacun), par-dessus le Mesh, avec la couleur RGBA de chaque
    /// sommet (l'alpha mélange au Mesh) ; vide : aucun.
    pub fn set_highlight(&mut self, triangles: &[[f32; 3]], colors: &[[u8; 4]]) {
        self.highlight = (!triangles.is_empty()).then(|| {
            let buffer = |contents: &[u8]| {
                self.device
                    .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: None,
                        contents,
                        usage: wgpu::BufferUsages::VERTEX,
                    })
            };
            (
                buffer(bytemuck::cast_slice(triangles)),
                buffer(bytemuck::cast_slice(colors)),
                triangles.len() as u32,
            )
        });
        self.request_redraw();
    }

    /// (min, max) du champ, qui recevront les deux bouts de la colormap.
    pub fn set_range(&mut self, range: (f64, f64)) {
        self.range = range;
        self.write_uniforms();
    }

    pub fn set_theme(&mut self, theme: Theme) {
        if self.theme != theme {
            self.theme = theme;
            self.write_uniforms();
        }
    }

    pub fn set_camera(&mut self, camera: &Camera) {
        self.camera = camera.matrix(self.size());
        let (_, _, eye) = camera.axes();
        self.eye = eye.map(|c| c as f32);
        self.write_uniforms();
    }

    fn write_uniforms(&mut self) {
        let mut data = [0.0_f32; 28];
        data[..16].copy_from_slice(self.camera.as_flattened());
        data[16] = self.range.0 as f32;
        data[17] = self.range.1 as f32;
        data[20..23].copy_from_slice(&self.eye);
        data[24..28].copy_from_slice(&self.theme.edge());
        self.queue
            .write_buffer(&self.uniforms, 0, bytemuck::cast_slice(&data));
        self.request_redraw();
    }

    pub fn resize(&mut self, width: u32, height: u32) {
        // fenêtre réduite : taille nulle, refusée par wgpu
        if width == 0 || height == 0 {
            return;
        }
        match &mut self.target {
            Target::Window {
                surface, config, ..
            } => {
                config.width = width;
                config.height = height;
                surface.configure(&self.device, config);
            }
            Target::Offscreen {
                texture,
                width: w,
                height: h,
                ..
            } => {
                if (*w, *h) == (width, height) {
                    return;
                }
                *texture = offscreen_texture(&self.device, width, height);
                (*w, *h) = (width, height);
            }
        }
        self.depth = depth_view(&self.device, width, height);
        // la matrice dépend du rapport de la fenêtre : l'appelant renvoie la caméra
        self.write_uniforms();
    }

    /// Facteur d'échelle de l'écran d'une texture partagée (écran Retina : 2) ; sans effet sur une
    /// fenêtre, qui le lit d'elle-même.
    pub fn set_scale_factor(&mut self, factor: f64) {
        if let Target::Offscreen { scale, .. } = &mut self.target {
            *scale = factor;
        }
    }

    /// La texture d'une cible hors écran (format sRGB 8 bits), que `draw` remplit ; `None` pour une
    /// fenêtre. Elle change quand `resize` change la taille.
    pub fn texture(&self) -> Option<&wgpu::Texture> {
        match &self.target {
            Target::Offscreen { texture, .. } => Some(texture),
            Target::Window { .. } => None,
        }
    }

    /// Dessine dans la texture d'une cible hors écran, sans relire les pixels.
    pub fn draw(&mut self) {
        let Target::Offscreen { texture, .. } = &self.target else {
            return;
        };
        let view = texture.create_view(&Default::default());
        let encoder = self.encode(&view);
        self.queue.submit([encoder.finish()]);
        self.text.atlas.trim();
    }

    /// Dessine dans la fenêtre.
    pub fn render(&mut self) {
        let Target::Window {
            window,
            surface,
            config,
        } = &self.target
        else {
            return;
        };
        let frame = match surface.get_current_texture() {
            wgpu::CurrentSurfaceTexture::Success(frame)
            | wgpu::CurrentSurfaceTexture::Suboptimal(frame) => frame,
            wgpu::CurrentSurfaceTexture::Outdated | wgpu::CurrentSurfaceTexture::Lost => {
                surface.configure(&self.device, config);
                window.request_redraw();
                return;
            }
            // fenêtre masquée, délai dépassé… : image sautée
            _ => return,
        };
        let window = window.clone();
        let view = frame.texture.create_view(&Default::default());
        let encoder = self.encode(&view);
        self.queue.submit([encoder.finish()]);
        window.pre_present_notify();
        self.queue.present(frame);
        self.text.atlas.trim();
    }

    /// Dessine hors écran et rend les pixels, RGBA 8 bits, ligne par ligne du haut vers le bas.
    pub fn render_pixels(&mut self) -> Result<Vec<u8>, Box<dyn Error>> {
        let Target::Offscreen {
            texture,
            width,
            height,
            ..
        } = &self.target
        else {
            return Err("capture demandée sur une fenêtre".into());
        };
        let (texture, width, height) = (texture.clone(), *width, *height);
        let view = texture.create_view(&Default::default());
        let mut encoder = self.encode(&view);
        // les lignes de la copie sont alignées sur 256 octets
        let padded = (width * 4).div_ceil(wgpu::COPY_BYTES_PER_ROW_ALIGNMENT)
            * wgpu::COPY_BYTES_PER_ROW_ALIGNMENT;
        let readback = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size: u64::from(padded * height),
            usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
            mapped_at_creation: false,
        });
        encoder.copy_texture_to_buffer(
            texture.as_image_copy(),
            wgpu::TexelCopyBufferInfo {
                buffer: &readback,
                layout: wgpu::TexelCopyBufferLayout {
                    offset: 0,
                    bytes_per_row: Some(padded),
                    rows_per_image: Some(height),
                },
            },
            wgpu::Extent3d {
                width,
                height,
                depth_or_array_layers: 1,
            },
        );
        self.queue.submit([encoder.finish()]);
        let (sender, receiver) = std::sync::mpsc::channel();
        readback
            .slice(..)
            .map_async(wgpu::MapMode::Read, move |result| {
                let _ = sender.send(result);
            });
        self.device.poll(wgpu::PollType::wait_indefinitely())?;
        receiver.recv()??;
        let mapped = readback.slice(..).get_mapped_range()?;
        let row = (width * 4) as usize;
        let mut pixels = Vec::with_capacity(row * height as usize);
        for line in mapped.chunks(padded as usize).take(height as usize) {
            pixels.extend_from_slice(&line[..row]);
        }
        Ok(pixels)
    }

    /// Enregistre le dessin de l'image dans `view` : Mesh, contour, légende, texte.
    fn encode(&mut self, view: &wgpu::TextureView) -> wgpu::CommandEncoder {
        let text_ready = self.prepare_text();
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(self.theme.background()),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &self.depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Discard,
                    }),
                    stencil_ops: None,
                }),
                ..Default::default()
            });
            pass.set_bind_group(0, &self.bind_group, &[]);
            if let Some(mesh) = &self.mesh {
                // remplissage d'abord, arêtes par-dessus
                pass.set_pipeline(&self.fill_pipeline);
                pass.set_vertex_buffer(0, mesh.positions.slice(..));
                pass.set_vertex_buffer(1, mesh.values.slice(..));
                pass.set_index_buffer(mesh.triangles.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..mesh.triangle_indices, 0, 0..1);
                if let Some((positions, colors, count)) = &self.highlight {
                    pass.set_pipeline(&self.highlight_pipeline);
                    pass.set_vertex_buffer(0, positions.slice(..));
                    pass.set_vertex_buffer(1, colors.slice(..));
                    pass.draw(0..*count, 0..1);
                }
                pass.set_pipeline(&self.edge_pipeline);
                pass.set_vertex_buffer(0, mesh.lines.slice(..));
                pass.draw(0..mesh.line_vertices, 0..1);
            }
            // légende par-dessus le Mesh
            if let Some(legend) = &self.legend {
                pass.set_vertex_buffer(0, legend.positions.slice(..));
                pass.set_vertex_buffer(1, legend.values.slice(..));
                pass.set_pipeline(&self.legend_fill_pipeline);
                pass.set_index_buffer(legend.bar.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..6, 0, 0..1);
                pass.set_pipeline(&self.legend_line_pipeline);
                pass.set_index_buffer(legend.lines.slice(..), wgpu::IndexFormat::Uint32);
                pass.draw_indexed(0..legend.line_indices, 0, 0..1);
            }
            if text_ready {
                let t = &self.text;
                if let Err(e) = t.renderer.render(&t.atlas, &t.viewport, &mut pass) {
                    eprintln!("légende : texte non affiché ({e})");
                }
            }
        }
        encoder
    }

    /// Envoie les glyphes de la légende au GPU ; faux (texte sauté pour cette image) en cas d'échec.
    fn prepare_text(&mut self) -> bool {
        let size = self.size();
        let t = &mut self.text;
        t.viewport.update(
            &self.queue,
            glyphon::Resolution {
                width: size[0] as u32,
                height: size[1] as u32,
            },
        );
        let bounds = glyphon::TextBounds {
            left: 0,
            top: 0,
            right: size[0] as i32,
            bottom: size[1] as i32,
        };
        let areas = t.blocks.iter().map(|(buffer, at)| glyphon::TextArea {
            buffer,
            left: at[0],
            top: at[1],
            scale: 1.0,
            bounds,
            default_color: self.theme.text(),
            custom_glyphs: &[],
        });
        match t.renderer.prepare(
            &self.device,
            &self.queue,
            &mut t.fonts,
            &mut t.atlas,
            &t.viewport,
            areas,
            &mut t.swash_cache,
        ) {
            Ok(()) => true,
            Err(e) => {
                eprintln!("légende : texte non affiché ({e})");
                false
            }
        }
    }
}

/// Bind group du shader : uniforms, LUT de `colormap` (texture 1D sRGB, lue sans interpolation pour
/// des bandes nettes) et son sampler.
fn lut_bind_group(
    device: &wgpu::Device,
    queue: &wgpu::Queue,
    layout: &wgpu::BindGroupLayout,
    uniforms: &wgpu::Buffer,
    sampler: &wgpu::Sampler,
    colormap: Colormap,
) -> wgpu::BindGroup {
    let lut = colormap.lut();
    let texture = device.create_texture_with_data(
        queue,
        &wgpu::TextureDescriptor {
            label: None,
            size: wgpu::Extent3d {
                width: lut.len() as u32,
                height: 1,
                depth_or_array_layers: 1,
            },
            mip_level_count: 1,
            sample_count: 1,
            dimension: wgpu::TextureDimension::D1,
            format: wgpu::TextureFormat::Rgba8UnormSrgb,
            usage: wgpu::TextureUsages::TEXTURE_BINDING,
            view_formats: &[],
        },
        wgpu::util::TextureDataOrder::LayerMajor,
        lut.as_flattened(),
    );
    device.create_bind_group(&wgpu::BindGroupDescriptor {
        label: None,
        layout,
        entries: &[
            wgpu::BindGroupEntry {
                binding: 0,
                resource: uniforms.as_entire_binding(),
            },
            wgpu::BindGroupEntry {
                binding: 1,
                resource: wgpu::BindingResource::TextureView(
                    &texture.create_view(&Default::default()),
                ),
            },
            wgpu::BindGroupEntry {
                binding: 2,
                resource: wgpu::BindingResource::Sampler(sampler),
            },
        ],
    })
}
