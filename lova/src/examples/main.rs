use std::io::Result;
use lattirust_arithmetic::ring::Z2_64;
use lova::util::{OptimizationMode, PublicParameters};

mod tui {
    use std::io::{self, stdout, Stdout};

    use crossterm::event::{Event, KeyCode, KeyEvent, KeyEventKind};
    use crossterm::{event, execute, terminal::*};
    use log::LevelFilter;
    use ratatui::prelude::*;
    use ratatui::symbols::border;
    use ratatui::widgets::block::Position;
    use ratatui::widgets::block::Title;
    use ratatui::widgets::Borders;
    use ratatui::widgets::Widget;
    use ratatui::widgets::{Block, Padding, Paragraph};
    use tui_logger::TuiLoggerWidget;
    use tui_logger::*;
    use tui_logger::{TuiLoggerSmartWidget, TuiWidgetState};

    use lova::util::header;

    use crate::tui;

    use super::*;

    /// A type alias for the terminal type used in this application
    pub type Tui = Terminal<CrosstermBackend<Stdout>>;

    /// Initialize the terminal
    pub fn init() -> io::Result<Tui> {
        // Init tui_logger
        tui_logger::init_logger(log::LevelFilter::Debug).unwrap();
        tui_logger::set_default_level(log::LevelFilter::Debug);

        // Init terminal
        execute!(stdout(), EnterAlternateScreen)?;
        enable_raw_mode()?;
        Terminal::new(CrosstermBackend::new(stdout()))
    }

    /// Restore the terminal to its original state
    pub fn restore() -> io::Result<()> {
        execute!(stdout(), LeaveAlternateScreen)?;
        disable_raw_mode()?;
        Ok(())
    }

    #[derive(Debug, Default)]
    pub struct App {
        exit: bool,
    }

    impl App {
        /// runs the application's main loop until the user quits
        pub fn run(&mut self, terminal: &mut tui::Tui) -> io::Result<()> {
            while !self.exit {
                terminal.draw(|frame| self.render_frame(frame))?;
                self.handle_events()?;
            }
            Ok(())
        }

        fn render_frame(&self, frame: &mut Frame) {
            let header = header();
            let header_width = header
                .lines()
                .map(|line| line.chars().count())
                .max()
                .unwrap_or(0) as u16;
            let header_height = header.lines().collect::<Vec<_>>().len() as u16;

            let outer_layout = Layout::default()
                .direction(Direction::Vertical)
                .constraints(vec![Constraint::Length(header_height), Constraint::Min(1)])
                .split(frame.size());

            let inner_layout = Layout::default()
                .direction(Direction::Horizontal)
                .constraints(vec![Constraint::Length(header_width), Constraint::Min(1)])
                .split(outer_layout[0]);

            frame.render_widget(
                Paragraph::new(header)
                    .centered()
                    .block(Block::new().style(Style::new()).padding(Padding::new(
                        0,
                        0,
                        (inner_layout[0].height - header_height) / 2,
                        0,
                    ))),
                inner_layout[0],
            );
            frame.render_widget(
                Paragraph::new("htop").block(Block::new().borders(Borders::ALL)),
                inner_layout[1],
            );

            let filter_state = TuiWidgetState::new()
                .set_default_display_level(LevelFilter::Debug);
            frame.render_widget(
                TuiLoggerSmartWidget::default().style_error(Style::default().fg(Color::Red))
                    .style_debug(Style::default().fg(Color::Green))
                    .style_warn(Style::default().fg(Color::Yellow))
                    .style_trace(Style::default().fg(Color::Magenta))
                    .style_info(Style::default().fg(Color::Cyan))
                    .output_separator(':')
                    .output_timestamp(Some("%H:%M:%S".to_string()))
                    .output_level(Some(TuiLoggerLevelOutput::Abbreviated))
                    .output_target(true)
                    .output_file(true)
                    .output_line(true)
                    .state(&filter_state),
                outer_layout[1],
            );
        }

        fn handle_events(&mut self) -> io::Result<()> {
            match event::read()? {
                // it's important to check that the event is a key press event as
                // crossterm also emits key release and repeat events on Windows.
                Event::Key(key_event) if key_event.kind == KeyEventKind::Press => {
                    self.handle_key_event(key_event)
                }
                _ => {}
            };
            Ok(())
        }

        fn handle_key_event(&mut self, key_event: KeyEvent) {
            match key_event.code {
                KeyCode::Char('q') => self.exit(),
                _ => {}
            }
        }

        fn exit(&mut self) {
            self.exit = true;
        }
    }

    impl Widget for &App {
        fn render(self, area: Rect, buf: &mut Buffer) {
            let title = Title::from(" Counter App Tutorial ".bold());
            let instructions = Title::from(Line::from(vec![" Quit ".into(), "<Q> ".blue().bold()]));
            let block = Block::default()
                .title(title.alignment(Alignment::Center))
                .title(
                    instructions
                        .alignment(Alignment::Center)
                        .position(Position::Bottom),
                )
                .borders(Borders::ALL)
                .border_set(border::THICK);

            let counter_text = Text::from(vec![Line::from(vec!["Value: ".into(), "...".yellow()])]);

            Paragraph::new(counter_text)
                .centered()
                .block(block)
                .render(area, buf);
        }
    }
}
type F = Z2_64;

const N: usize = 1 << 16;

fn main() -> Result<()> {
    use tui::*;
    let mut terminal = tui::init()?;
    let app_result = App::default().run(&mut terminal);
    tui::restore()?;
    app_result
}
