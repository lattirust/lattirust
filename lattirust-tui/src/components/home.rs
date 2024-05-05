use std::sync::{Arc, Mutex, RwLock};
use std::time::Duration;

use color_eyre::eyre::Result;
use ratatui::{prelude::*, widgets::*};
use tokio::sync::mpsc::UnboundedSender;

use crate::{action::Action, config::Config};

use super::{Component, Frame};

pub struct Home<Instance, Witness, PublicParameters, Proof> {
    command_tx: Option<UnboundedSender<Action>>,
    config: Config,
    instance: Arc<RwLock<Option<Instance>>>,
    witness: Arc<RwLock<Option<Witness>>>,
    public_parameters: Arc<RwLock<Option<PublicParameters>>>,
    proof: Arc<RwLock<Option<Proof>>>,
}

impl<Instance, Witness, PublicParameters, Proof> Default
    for Home<Instance, Witness, PublicParameters, Proof>
{
    fn default() -> Self {
        Self {
            command_tx: None,
            config: Config::default(),
            instance: Arc::new(RwLock::new(None)),
            witness: Arc::new(RwLock::new(None)),
            public_parameters: Arc::new(RwLock::new(None)),
            proof: Arc::new(RwLock::new(None)),
        }
    }
}

impl<Instance, Witness, Proof> Home<Instance, Witness, usize, Proof> {
    pub fn new() -> Self {
        Self::default()
    }

    async fn setup(public_parameters: Arc<RwLock<Option<usize>>>, tx: UnboundedSender<Action>) {
        tokio::time::sleep(Duration::from_secs(2)).await;
        if let Ok(mut guard) = public_parameters.try_write() {
            *guard = match *guard {
                None => Some(0).into(),
                Some(x) => Some(x + 1).into(),
            };
        } else {
            println!("Couldn't get write access, sorry!")
        };

        tx.send(Action::Prove).unwrap();
    }
}

impl<Instance: 'static, Witness: 'static, Proof: 'static> Component
    for Home<Instance, Witness, usize, Proof>
{
    fn register_action_handler(&mut self, tx: UnboundedSender<Action>) -> Result<()> {
        self.command_tx = Some(tx);
        Ok(())
    }

    fn register_config_handler(&mut self, config: Config) -> Result<()> {
        self.config = config;
        Ok(())
    }

    fn update(&mut self, action: Action) -> Result<Option<Action>> {
        match action {
            Action::Tick => {}
            Action::Setup => {
                assert!(self.command_tx.is_some());
                let tx = self.command_tx.clone().unwrap();
                let c = self.public_parameters.clone();
                tokio::spawn(Self::setup(c, tx));
            }
            _ => {}
        }
        Ok(None)
    }

    fn draw(&mut self, f: &mut Frame<'_>, area: Rect) -> Result<()> {
        f.render_widget(
            Paragraph::new(format!("hello world: {:?}", *self.public_parameters)),
            area,
        );
        Ok(())
    }
}
