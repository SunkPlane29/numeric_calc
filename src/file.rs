use std::fs::File;
use std::io::prelude::*;

//TODO: find a way to export this
pub struct DataFile {
    f: File,
}

impl DataFile {
    pub fn create(path: &str) -> DataFile {
        let file = File::options()
            .write(true)
            .append(false)
            .create(true)
            .open(path).unwrap();

        DataFile{f: file}
    }

    //TODO: create more write methods, and possibly alse read methods
    pub fn  write(&mut self, first_column: f64, second_column: f64) {
        let row = format!(" {} {}\n", first_column, second_column);
        self.f.write_all(row.as_bytes()).unwrap();
    }
}
