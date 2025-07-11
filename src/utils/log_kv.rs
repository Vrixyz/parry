#[macro_export]
/// log only if `log_kv_serde` feature is enabled, in order to be able to pass key values.
/// First argument is the debug level, then regular `log` arguments.
///
/// Example usage:
///
/// ```rs
/// log_kv!(debug, key1:serde = 1, key2:serde = "a string"; "Some description");
/// ```
macro_rules! log_kv {
    ($level:ident, $($arg:tt)*) => {
        #[cfg(feature = "log_kv_serde")]
        {
            log::$level!($($arg)*);
        }

        #[cfg(not(feature = "log_kv_serde"))]
        {
            // Do nothing
        }
    };
}
