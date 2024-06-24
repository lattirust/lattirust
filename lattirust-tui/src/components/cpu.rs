use std::cell::Cell;
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::prelude::{Color, Constraint, Direction, Layout, Style, Text, Widget};
use ratatui::widgets::{Block, Borders, Gauge, Paragraph, Row, Table};
use std::sync::{Arc, Mutex};
use color_eyre::owo_colors::OwoColorize;
use sysinfo::System;

#[derive(Clone)]
pub struct CPUWidget {
    system: Arc<Mutex<System>>,
}

impl Default for CPUWidget {
    fn default() -> Self {
        CPUWidget {
            system: Arc::new(Mutex::new(System::new_all())),
        }
    }
}

impl Widget for CPUWidget {
    fn render(mut self, area: Rect, buf: &mut Buffer)
    where
        Self: Sized,
    {
        let mut system = self.system.lock().unwrap();
        system.refresh_cpu();

        let cpus = system.cpus();

        let avg_cpu_utilization = cpus.iter().map(|cpu| cpu.cpu_usage()).sum::<f32>() / cpus.len() as f32;
        let cpu_color = if avg_cpu_utilization < 100.0 / cpus.len() as f32 {
            Color::Green
        } else if avg_cpu_utilization < 50.0 {
            Color::Yellow
        } else {
            Color::Red
        };
        let cpu_gauge = Gauge::default()
            .gauge_style(Style::default().fg(cpu_color)).label(format!("CPU ({}%)", avg_cpu_utilization.round()))
            .percent(avg_cpu_utilization as u16);

        let used_memory = system.used_memory();
        let total_memory = system.total_memory();
        let memory_ratio = used_memory as f64 / total_memory as f64;
        let memory_gauge = Gauge::default()
            .gauge_style(Style::default().fg(Color::Blue))
            .label(format!("Memory ({}%)", memory_ratio.round()))
            .ratio(memory_ratio);

        let layout = Layout::default()
            .direction(Direction::Vertical)
            .constraints(
                [
                    Constraint::Percentage(50),
                    Constraint::Percentage(50),
                ]
                .as_ref(),
            )
            .margin(1)
            .split(area);
        let block = Block::default().borders(Borders::ALL).title("Perf");

        block.render(area, buf);
        cpu_gauge.render(layout[0], buf);
        memory_gauge.render(layout[1], buf);
    }
}
