use std::mem;

fn main() {
    println!("HI!")
    let w = 4;
    let hit_thres = .9;
    let delta = 10.;
    let hsp_thres = 0.;
    let S = ('A', 'T', 'G', 'C')

    let d_fasta: &'static str = "data/raw/chr22.maf.ancestors.42000000.complete.boreo.fa"
    let d_conf = get_d_conf()

    let d = get_prob_seq()
}

fn get_prob_seq(fasta: String, conf: &[f64], S: ) -> *{

}

fn seq_from_fasta(fasta: String, , )