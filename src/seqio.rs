pub fn fasta_reader(file: &str) -> (Vec<String>, Vec<String>) {
    let mut seqs: Vec<Vec<String>> = Vec::new();
    let mut headers: Vec<String> = Vec::new();

    let file = std::fs::read_to_string(file).expect("Failed to read file");
    for line in file.lines() {
        if line.starts_with('>') {
            headers.push(line.to_string());
            seqs.push(Vec::new());
        } else {
            seqs.last_mut().unwrap().push(line.to_string());
        }
    }

    let seqs_final = seqs
        .iter_mut()
        .map(|seq_list| seq_list.concat().to_uppercase())
        .collect::<Vec<String>>();

    (headers, seqs_final)
}

pub fn remove_end_insertions(mut seq_array: Vec<Vec<u8>>) -> Vec<Vec<u8>> {
    for seq in seq_array.iter_mut() {
        // Remove the right ends
        for base in seq.iter_mut() {
            match base {
                b'-' => {
                    *base = b' ';
                }
                _ => break,
            }
        }
        // Remove the left ends
        for base in seq.iter_mut().rev() {
            match base {
                b'-' => {
                    *base = b' ';
                }
                _ => break,
            }
        }
    }
    seq_array
}
