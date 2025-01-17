#[derive(Debug)]
pub enum InsecureOptionError {
    /// Field Size is not big enough
    FieldSize,
    /// Number of security bits is not enough
    LowSecurityBits,
}
