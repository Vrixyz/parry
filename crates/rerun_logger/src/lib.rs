use std::sync::atomic::AtomicU32;

use log::{Log, Record};
pub mod parry;

pub struct RerunLogger;

pub static RERUN_LOGGER: RerunLogger = RerunLogger;

static RR_CELL: std::sync::OnceLock<rerun::RecordingStream> = std::sync::OnceLock::new();
static RECORD: std::sync::OnceLock<RecordFn> = std::sync::OnceLock::new();
pub static TIME: AtomicU32 = AtomicU32::new(0);

pub fn init_rr() {
    let _ = RR_CELL.get_or_init(|| {
        rerun::RecordingStreamBuilder::new("rerun_test_311_convex_hull")
            .spawn()
            .unwrap()
    });

    get_rr().set_duration_secs("sim_time", 0);
}

pub fn get_rr() -> &'static rerun::RecordingStream {
    RR_CELL.get().unwrap()
}
pub type RecordFn = Option<Box<dyn Fn(&Record) + Sync + Send>>;

pub fn init_record(record: RecordFn) {
    RECORD.get_or_init(|| record);
}
pub fn get_record() -> &'static RecordFn {
    RECORD.get_or_init(|| None)
}

impl Log for RerunLogger {
    fn enabled(&self, _metadata: &log::Metadata) -> bool {
        true
    }

    fn log(&self, record: &log::Record) {
        if !self.enabled(record.metadata()) {
            return;
        }
        std::println!(
            "{}:{} -- {}",
            record.level(),
            record.target(),
            record.args()
        );

        if let Some(handle_record) = &get_record() {
            handle_record(record);
        }
        //self.forward_logger.log(record);
    }
    fn flush(&self) {}
}
