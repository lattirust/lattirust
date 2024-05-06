use std::sync::{Arc, RwLock};
use std::time::Duration;

use color_eyre::eyre::Result;
use color_eyre::owo_colors::XtermColors;
use color_eyre::owo_colors::XtermColors::ElectricPurple;
use log::{info, LevelFilter};
use ratatui::{prelude::*, widgets::*};
use tokio::sync::mpsc::UnboundedSender;
use tracing::instrument::WithSubscriber;
use tui_logger::{TuiLoggerLevelOutput, TuiLoggerSmartWidget, TuiLoggerWidget, TuiWidgetState};

use crate::{action::Action, config::Config};

use super::{Component, Frame};

const THEME_COLOR: Color = Color::Magenta;

pub struct Home<Instance, Witness, PublicParameters, Proof> {
    command_tx: Option<UnboundedSender<Action>>,
    config: Config,
    instance: Arc<RwLock<Option<Instance>>>,
    witness: Arc<RwLock<Option<Witness>>>,
    public_parameters: Arc<RwLock<Option<PublicParameters>>>,
    proof: Arc<RwLock<Option<Proof>>>,
    header: String,
    header_widget: Paragraph<'static>,
    header_layout: Rect,
    perf_layout: Rect,
    logs_layout: Rect,
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
            header: String::default(),
            header_widget: Paragraph::new(String::default()),
            header_layout: Rect::default(),
            perf_layout: Rect::default(),
            logs_layout: Rect::default(),
        }
    }
}

impl<Instance, Witness, Proof> Home<Instance, Witness, usize, Proof> {
    pub fn new(header: String) -> Self {
        let mut res = Self::default();
        res.header = header;
        res.set_layout(Rect::default());
        res
    }

    fn set_layout(&mut self, area: Rect) {
        let header_width = self
            .header
            .lines()
            .map(|line| line.chars().count())
            .max()
            .unwrap_or(0) as u16;
        let header_height = self.header.lines().collect::<Vec<_>>().len() as u16;

        let outer_layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints(vec![Constraint::Length(header_height), Constraint::Min(1)])
            .split(area);

        let inner_layout = Layout::default()
            .direction(Direction::Horizontal)
            .constraints(vec![Constraint::Length(header_width), Constraint::Min(1)])
            .split(outer_layout[0]);
        self.logs_layout = outer_layout[1];
        self.header_layout = inner_layout[0];
        self.perf_layout = inner_layout[1];

        let padding_top = if self.header_layout.height < header_height {
            0
        } else {
            (self.header_layout.height - header_height) / 2
        };
        self.header_widget = Paragraph::new(self.header.clone()).centered().block(
            Block::new()
                .style(Style::new())
                .padding(Padding::top(padding_top)),
        );
    }

    #[tracing::instrument]
    fn test() {
        std::thread::sleep(Duration::from_secs(1));
        tracing::debug!("Test");
    }

    #[tracing::instrument]
    async fn setup(public_parameters: Arc<RwLock<Option<usize>>>, tx: UnboundedSender<Action>) {
        tokio::time::sleep(Duration::from_secs(2)).await;
        if let Ok(mut guard) = public_parameters.try_write() {
            info!("[log] Setting up public parameters");
            tracing::event!(
                tracing::Level::DEBUG,
                "[tracing::event] Setting up public parameters"
            );
            tracing::trace!("[tracing] Setting up public parameters");
            *guard = match *guard {
                None => Some(0).into(),
                Some(x) => Some(x + 1).into(),
            };
            Self::test();
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

    fn draw(&mut self, frame: &mut Frame<'_>, area: Rect) -> Result<()> {
        self.set_layout(frame.size());

        frame.render_widget(self.header_widget.clone(), self.header_layout);
        frame.render_widget(
            Paragraph::new(format!("htop:  {:?}", self.public_parameters))
                .block(Block::new().borders(Borders::ALL)),
            self.perf_layout,
        );

        let filter_state = TuiWidgetState::new().set_default_display_level(LevelFilter::Trace);
        frame.render_widget(
            TuiLoggerWidget::default()
                .style_error(Style::default().fg(Color::Red))
                .style_warn(Style::default().fg(Color::Yellow))
                .style_info(Style::default().fg(Color::Green))
                .style_debug(Style::default().fg(Color::Magenta))
                .style_trace(Style::default().fg(Color::Gray))
                .output_separator(' ')
                .output_timestamp(Some("%H:%M:%S".to_string()))
                .output_level(Some(TuiLoggerLevelOutput::Long))
                .output_target(false)
                .output_file(false)
                .output_line(false)
                .state(&filter_state)
                .block(
                    Block::new()
                        .borders(Borders::ALL)
                        .title_top(" Logs ")
                        .title_style(Style::default().fg(THEME_COLOR)),
                ),
            self.logs_layout,
        );
        Ok(())
    }
}
