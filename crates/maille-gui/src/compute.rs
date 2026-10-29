//! Le calcul en tâche de fond : `maille_io::solving::solve_file` (le calcul de `maille solve`) tourne
//! dans un thread, l'interface suit sa progression par un canal et peut l'annuler. L'état de la
//! progression (`ComputeState`) est pur et testé ; `Job` ne fait que les fils et le thread.

use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicBool, Ordering};
use std::sync::mpsc::{self, Receiver, TryRecvError};
use std::sync::Arc;
use std::thread::{self, JoinHandle};

use maille_core::model::{Event, Flow};
use maille_core::MailleError;

/// Ce que le thread dit à l'interface.
#[derive(Debug)]
pub enum Message {
    /// Une étape du calcul (maillage fabriqué, calcul thermique, calcul…).
    Log(String),
    Progress(Event),
    /// Le calcul est fini : les fichiers écrits (`.pvd` en dernier), l'arrêt demandé, ou l'échec.
    Done(Outcome),
}

#[derive(Debug, PartialEq)]
pub enum Outcome {
    Finished(Vec<PathBuf>),
    Cancelled,
    Failed(String),
}

/// Où en est le calcul, tel que l'interface l'affiche.
#[derive(Debug, Clone, Default, PartialEq)]
pub enum ComputeState {
    #[default]
    Idle,
    Running(Running),
    /// Fini : le `.pvd` des résultats.
    Finished(PathBuf),
    Cancelled,
    Failed(String),
}

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Running {
    /// La dernière étape annoncée.
    pub phase: String,
    /// Increments finis et total de la Step (0 sur 0 : calcul sans Step, ou pas encore commencé).
    pub done: usize,
    pub total: usize,
    /// Dernière itération de Newton-Raphson et son résidu relatif.
    pub iteration: Option<(usize, f64)>,
    /// Coupes d'incrément faites.
    pub cuts: usize,
    /// L'arrêt est demandé : le calcul s'arrête à la prochaine itération.
    pub stopping: bool,
}

impl ComputeState {
    pub fn is_running(&self) -> bool {
        matches!(self, ComputeState::Running(_))
    }

    /// Prend en compte un message du thread.
    pub fn apply(&mut self, message: Message) {
        match message {
            Message::Log(phase) => {
                if let ComputeState::Running(r) = self {
                    r.phase = phase;
                    r.done = 0;
                    r.total = 0;
                    r.iteration = None;
                }
            }
            Message::Progress(event) => {
                if let ComputeState::Running(r) = self {
                    match event {
                        Event::Increment { done, total, .. } => {
                            r.done = done;
                            r.total = total;
                        }
                        Event::Iteration {
                            iteration,
                            residual,
                            ..
                        } => r.iteration = Some((iteration, residual)),
                        Event::Cut { cuts, .. } => r.cuts = r.cuts.max(cuts),
                    }
                }
            }
            Message::Done(outcome) => {
                *self = match outcome {
                    Outcome::Finished(files) => match files.into_iter().last() {
                        Some(pvd) => ComputeState::Finished(pvd),
                        None => ComputeState::Failed("le calcul n'a rien écrit".into()),
                    },
                    Outcome::Cancelled => ComputeState::Cancelled,
                    Outcome::Failed(why) => ComputeState::Failed(why),
                };
            }
        }
    }

    /// Part du calcul faite, de 0 à 1, quand on la connaît.
    pub fn fraction(&self) -> Option<f32> {
        match self {
            ComputeState::Running(r) if r.total > 0 => Some(r.done as f32 / r.total as f32),
            _ => None,
        }
    }
}

/// Un calcul lancé dans un thread.
pub struct Job {
    receiver: Receiver<Message>,
    cancel: Arc<AtomicBool>,
    handle: Option<JoinHandle<()>>,
    finished: bool,
}

impl Job {
    /// Lance `maille solve` sur ce modèle. `wake` est appelé à chaque message, pour que l'interface se
    /// redessine.
    pub fn start(model: &Path, wake: impl Fn() + Send + 'static) -> Self {
        let (sender, receiver) = mpsc::channel();
        let cancel = Arc::new(AtomicBool::new(false));
        let flag = Arc::clone(&cancel);
        let model = model.to_path_buf();
        let handle = thread::spawn(move || {
            let send = |message: Message| {
                // l'interface a pu fermer : plus personne n'écoute
                let _ = sender.send(message);
                wake();
            };
            let result = maille_io::solving::solve_file(
                &model,
                None,
                &mut |event| {
                    send(Message::Progress(event.clone()));
                    if flag.load(Ordering::Relaxed) {
                        Flow::Stop
                    } else {
                        Flow::Continue
                    }
                },
                &mut |line| send(Message::Log(line.to_string())),
            );
            send(Message::Done(match result {
                Ok(files) => Outcome::Finished(files),
                Err(e) if matches!(e.downcast_ref(), Some(MailleError::Cancelled)) => {
                    Outcome::Cancelled
                }
                Err(e) => Outcome::Failed(e.to_string()),
            }));
        });
        Self {
            receiver,
            cancel,
            handle: Some(handle),
            finished: false,
        }
    }

    /// Demande l'arrêt : le calcul s'arrête à son prochain rappel.
    pub fn cancel(&self) {
        self.cancel.store(true, Ordering::Relaxed);
    }

    /// Les messages arrivés depuis le dernier appel ; le dernier est `Done`.
    pub fn poll(&mut self) -> Vec<Message> {
        let mut out = Vec::new();
        loop {
            match self.receiver.try_recv() {
                Ok(message) => {
                    self.finished |= matches!(message, Message::Done(_));
                    out.push(message);
                }
                Err(TryRecvError::Empty) => break,
                Err(TryRecvError::Disconnected) => {
                    if !self.finished {
                        self.finished = true;
                        out.push(Message::Done(Outcome::Failed(
                            "le calcul s'est arrêté sans rendre de résultat".into(),
                        )));
                    }
                    break;
                }
            }
        }
        if self.finished {
            if let Some(handle) = self.handle.take() {
                let _ = handle.join();
            }
        }
        out
    }

    pub fn is_finished(&self) -> bool {
        self.finished
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{Duration, Instant};

    /// Copie ces exemples dans un dossier neuf (le calcul écrit à côté du modèle).
    fn sandbox(name: &str, files: &[&str]) -> PathBuf {
        let dir =
            std::env::temp_dir().join(format!("maille-compute-{name}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        for f in files {
            let from = Path::new(env!("CARGO_MANIFEST_DIR"))
                .join("../../examples")
                .join(f);
            fs::copy(from, dir.join(f)).unwrap();
        }
        dir
    }

    /// Attend la fin du calcul en appliquant ses messages ; rend l'état final et tous les messages.
    fn finish(job: &mut Job) -> (ComputeState, usize) {
        let mut state = ComputeState::Running(Running::default());
        let mut count = 0;
        let start = Instant::now();
        while !job.is_finished() {
            assert!(
                start.elapsed() < Duration::from_secs(120),
                "calcul trop long"
            );
            for m in job.poll() {
                count += 1;
                state.apply(m);
            }
            thread::sleep(Duration::from_millis(5));
        }
        (state, count)
    }

    #[test]
    fn a_job_solves_a_model_file_and_reports_the_pvd() {
        let dir = sandbox("traction", &["traction.toml", "plaque.msh"]);
        let mut job = Job::start(&dir.join("traction.toml"), || {});
        let (state, _) = finish(&mut job);
        let ComputeState::Finished(pvd) = state else {
            panic!("{state:?}");
        };
        assert_eq!(pvd, dir.join("traction.pvd"));
        assert!(pvd.exists() && dir.join("traction.vtu").exists());
    }

    #[test]
    fn a_job_gives_the_same_files_as_the_command_line_solver() {
        let a = sandbox("same-a", &["traction.toml", "plaque.msh"]);
        let b = sandbox("same-b", &["traction.toml", "plaque.msh"]);
        finish(&mut Job::start(&a.join("traction.toml"), || {}));
        maille_io::solving::solve_file(
            &b.join("traction.toml"),
            None,
            &mut |_| Flow::Continue,
            &mut |_| {},
        )
        .unwrap();
        for f in ["traction.vtu", "traction.pvd"] {
            assert_eq!(
                fs::read(a.join(f)).unwrap(),
                fs::read(b.join(f)).unwrap(),
                "{f}"
            );
        }
    }

    #[test]
    fn a_quasi_static_job_reports_its_progress_and_can_be_cancelled() {
        let dir = sandbox("hertz", &["hertz.toml", "hertz.msh"]);
        let mut job = Job::start(&dir.join("hertz.toml"), || {});
        job.cancel();
        let (state, messages) = finish(&mut job);
        assert_eq!(state, ComputeState::Cancelled);
        assert!(messages >= 2, "au moins une étape puis la fin");
        assert!(
            !dir.join("hertz.pvd").exists(),
            "rien n'est écrit après l'annulation"
        );
    }

    #[test]
    fn a_missing_model_fails_with_a_sentence() {
        let mut job = Job::start(Path::new("/nonexistent/model.toml"), || {});
        let (state, _) = finish(&mut job);
        let ComputeState::Failed(why) = state else {
            panic!("{state:?}");
        };
        assert!(why.contains("model.toml"), "{why}");
    }

    #[test]
    fn the_state_follows_the_messages() {
        let mut s = ComputeState::Running(Running::default());
        s.apply(Message::Log("calcul".into()));
        s.apply(Message::Progress(Event::Iteration {
            increment: 1,
            iteration: 2,
            residual: 1e-3,
        }));
        s.apply(Message::Progress(Event::Increment {
            done: 1,
            total: 4,
            time: 0.25,
        }));
        s.apply(Message::Progress(Event::Cut { time: 0.5, cuts: 1 }));
        let ComputeState::Running(r) = &s else {
            panic!()
        };
        assert_eq!((r.done, r.total, r.cuts), (1, 4, 1));
        assert_eq!(r.iteration, Some((2, 1e-3)));
        assert_eq!(s.fraction(), Some(0.25));
        s.apply(Message::Done(Outcome::Finished(vec![
            "a.vtu".into(),
            "a.pvd".into(),
        ])));
        assert_eq!(s, ComputeState::Finished("a.pvd".into()));
        assert!(!s.is_running());
    }
}
