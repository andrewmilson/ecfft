use core::fmt;
use std::error::Error;

#[derive(Debug, Copy, Clone, PartialEq, Eq)]
pub struct TryFromGeneralToShortWeierstrassCurveError;

impl fmt::Display for TryFromGeneralToShortWeierstrassCurveError {
    fn fmt(&self, fmt: &mut fmt::Formatter<'_>) -> fmt::Result {
        #[allow(deprecated)]
        self.description().fmt(fmt)
    }
}
impl Error for TryFromGeneralToShortWeierstrassCurveError {
    fn description(&self) -> &str {
        "conversion to short weierstrass without substitution attempted"
    }
}
