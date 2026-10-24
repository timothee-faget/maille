//! Contexte GPU (wgpu) attaché à la fenêtre : surface, pipelines, buffers du Mesh, dessin.

use std::error::Error;
use std::sync::Arc;

use wgpu::util::DeviceExt;
use winit::event_loop::ActiveEventLoop;
use winit::window::Window;

use crate::colormap::Colormap;
use crate::scene::Scene;
use crate::view::View;

/// Fond gris clair, sur lequel les couleurs de la colormap et le texte noir restent lisibles.
const BACKGROUND: wgpu::Color = wgpu::Color {
    r: 0.95,
    g: 0.95,
    b: 0.95,
    a: 1.0,
};

pub struct Renderer {
    window: Arc<Window>,
    surface: wgpu::Surface<'static>,
    device: wgpu::Device,
    queue: wgpu::Queue,
    config: wgpu::SurfaceConfiguration,
    fill_pipeline: wgpu::RenderPipeline,
    edge_pipeline: wgpu::RenderPipeline,
    bind_group_layout: wgpu::BindGroupLayout,
    sampler: wgpu::Sampler,
    /// Uniforms, LUT de la colormap et son sampler.
    bind_group: wgpu::BindGroup,
    /// Échelle, décalage et (min, max) du champ : voir `Uniforms` dans shader.wgsl.
    uniforms: wgpu::Buffer,
    /// Un [x, y] par Node.
    positions: wgpu::Buffer,
    /// Une valeur du champ par Node.
    values: wgpu::Buffer,
    triangles: wgpu::Buffer,
    triangle_indices: u32,
    edges: wgpu::Buffer,
    edge_indices: u32,
    view: View,
    range: (f64, f64),
}

impl Renderer {
    /// Prépare le GPU pour la géométrie de `scene` (triangles, arêtes) ; positions, valeurs et
    /// cadrage sont ensuite fixés par `set_positions`, `set_values` et `set_view`.
    pub fn new(
        event_loop: &ActiveEventLoop,
        window: Arc<Window>,
        scene: &Scene,
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
        let vertex_buffers = [
            Some(wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<[f32; 2]>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &wgpu::vertex_attr_array![0 => Float32x2],
            }),
            Some(wgpu::VertexBufferLayout {
                array_stride: std::mem::size_of::<f32>() as u64,
                step_mode: wgpu::VertexStepMode::Vertex,
                attributes: &wgpu::vertex_attr_array![1 => Float32],
            }),
        ];
        let pipeline = |topology, fragment| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: None,
                layout: Some(&layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some("vs_main"),
                    compilation_options: Default::default(),
                    buffers: &vertex_buffers,
                },
                primitive: wgpu::PrimitiveState {
                    topology,
                    ..Default::default()
                },
                depth_stencil: None,
                multisample: Default::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fragment),
                    compilation_options: Default::default(),
                    targets: &[Some(config.format.into())],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        let fill_pipeline = pipeline(wgpu::PrimitiveTopology::TriangleList, "fs_fill");
        let edge_pipeline = pipeline(wgpu::PrimitiveTopology::LineList, "fs_edge");

        let sampler = device.create_sampler(&wgpu::SamplerDescriptor::default());
        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            // 6 flottants utiles, complétés à 32 octets (multiple de 16 exigé par certains backends)
            size: std::mem::size_of::<[f32; 8]>() as u64,
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

        let node_count = scene.positions(0.0).len() as u64;
        let vertex_buffer = |bytes_per_node: u64| {
            device.create_buffer(&wgpu::BufferDescriptor {
                label: None,
                size: node_count * bytes_per_node,
                usage: wgpu::BufferUsages::VERTEX | wgpu::BufferUsages::COPY_DST,
                mapped_at_creation: false,
            })
        };
        let positions = vertex_buffer(std::mem::size_of::<[f32; 2]>() as u64);
        let values = vertex_buffer(std::mem::size_of::<f32>() as u64);
        let index_buffer = |indices: &[u32]| {
            device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
                label: None,
                contents: bytemuck::cast_slice(indices),
                usage: wgpu::BufferUsages::INDEX,
            })
        };
        let triangles = index_buffer(scene.triangles.as_flattened());
        let edges = index_buffer(scene.edges.as_flattened());

        Ok(Self {
            window,
            surface,
            device,
            queue,
            config,
            fill_pipeline,
            edge_pipeline,
            bind_group_layout,
            sampler,
            bind_group,
            uniforms,
            positions,
            values,
            triangle_indices: (scene.triangles.len() * 3) as u32,
            triangles,
            edge_indices: (scene.edges.len() * 2) as u32,
            edges,
            view: View::fit(scene.bounds(0.0), [1.0, 1.0]),
            range: (0.0, 1.0),
        })
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
        self.window.request_redraw();
    }

    pub fn set_title(&self, title: &str) {
        self.window.set_title(title);
    }

    /// Taille de la zone de dessin, en pixels.
    pub fn size(&self) -> [f64; 2] {
        [f64::from(self.config.width), f64::from(self.config.height)]
    }

    /// Positions des Nodes (une par Node, dans l'ordre de `Mesh::nodes`).
    pub fn set_positions(&mut self, positions: &[[f64; 2]]) {
        let data: Vec<[f32; 2]> = positions.iter().map(|p| p.map(|c| c as f32)).collect();
        self.queue
            .write_buffer(&self.positions, 0, bytemuck::cast_slice(&data));
        self.window.request_redraw();
    }

    /// Valeurs du champ aux Nodes, et (min, max) qui recevront les deux bouts de la colormap.
    pub fn set_values(&mut self, values: &[f64], range: (f64, f64)) {
        let data: Vec<f32> = values.iter().map(|&v| v as f32).collect();
        self.queue
            .write_buffer(&self.values, 0, bytemuck::cast_slice(&data));
        self.range = range;
        self.write_uniforms();
    }

    pub fn set_view(&mut self, view: View) {
        self.view = view;
        self.write_uniforms();
    }

    fn write_uniforms(&mut self) {
        let (scale, offset) = self.view.transform(self.size());
        let data = [
            scale[0],
            scale[1],
            offset[0],
            offset[1],
            self.range.0 as f32,
            self.range.1 as f32,
            0.0,
            0.0,
        ];
        self.queue
            .write_buffer(&self.uniforms, 0, bytemuck::cast_slice(&data));
        self.window.request_redraw();
    }

    pub fn resize(&mut self, width: u32, height: u32) {
        // fenêtre réduite : taille nulle, refusée par wgpu
        if width > 0 && height > 0 {
            self.config.width = width;
            self.config.height = height;
            self.surface.configure(&self.device, &self.config);
            self.write_uniforms();
        }
    }

    pub fn render(&mut self) {
        let frame = match self.surface.get_current_texture() {
            wgpu::CurrentSurfaceTexture::Success(frame)
            | wgpu::CurrentSurfaceTexture::Suboptimal(frame) => frame,
            wgpu::CurrentSurfaceTexture::Outdated | wgpu::CurrentSurfaceTexture::Lost => {
                self.surface.configure(&self.device, &self.config);
                self.window.request_redraw();
                return;
            }
            // fenêtre masquée, délai dépassé… : image sautée
            _ => return,
        };
        let view = frame.texture.create_view(&Default::default());
        let mut encoder = self.device.create_command_encoder(&Default::default());
        {
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &view,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(BACKGROUND),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                ..Default::default()
            });
            pass.set_bind_group(0, &self.bind_group, &[]);
            pass.set_vertex_buffer(0, self.positions.slice(..));
            pass.set_vertex_buffer(1, self.values.slice(..));
            // remplissage d'abord, arêtes par-dessus
            pass.set_pipeline(&self.fill_pipeline);
            pass.set_index_buffer(self.triangles.slice(..), wgpu::IndexFormat::Uint32);
            pass.draw_indexed(0..self.triangle_indices, 0, 0..1);
            pass.set_pipeline(&self.edge_pipeline);
            pass.set_index_buffer(self.edges.slice(..), wgpu::IndexFormat::Uint32);
            pass.draw_indexed(0..self.edge_indices, 0, 0..1);
        }
        self.queue.submit([encoder.finish()]);
        self.window.pre_present_notify();
        self.queue.present(frame);
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
