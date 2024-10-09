pub trait Transcript {
    fn append(&mut self, new_data: &[u8]);

    fn challenge(&mut self) -> [u8; 32];
}
