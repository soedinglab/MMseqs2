#[cfg(not(feature = "simd_avx2"))]
fn main() {}

#[cfg(feature = "simd_avx2")]
fn main() {
    use block_aligner::scan_block::*;
    use block_aligner::scores::*;
    use block_aligner::cigar::*;

    use image::{Rgb, RgbImage, ColorType};
    use image::codecs::png::{PngEncoder, CompressionType, FilterType};
    use imageproc::drawing::*;
    use imageproc::rect::Rect;

    use std::env;
    use std::io::BufWriter;
    use std::fs::File;

    let args = env::args().skip(1);

    let seqs = [
        // uc30_50_60
        (b"MVQATTWKKAIPGLSDEASSSPASELRAPLGGVRAMTMNELTRYSIKEPPSDELGSQLVNLYLQQLHTRYPFLDPAELWRLQKARTPVAHSESGNLSMTQRYGIFKLYMVFAIGATLLQLTNKSAEVSPERFYMTALQHMAAAKVPRTVQNIEAMTLLVVYHLRSASGLGLWYMIGLAMRTCIDLGLHRKNHERGLAPLVIQMHRRLFWTVYSLEIVIAISLGRPLSISERQIDVELPDTISVASVPCPSSPGETPVQPTSSNDNLQLANLLFQLRSIEARIHHSIYRTDKPLSALLPKLDKIYKQLEVWRLASIESLPPDGHVLDYPLLLYHRAVRMLIQPFMTILPVSDPYYVLCLRAAGSVCQMHKRLHQTIGYGHSFIAVQTIFVAGVTLLYGLWTQTHLVWSVTLADDLRACSLVLFVMSERAPWVRKYRDAFEVLVDAAMEKLRSGESSLAEMVAVAQTQAQAQSQSQGPRVGQFASGDETMRGPNPDTGPGSSSYGNGNGEHGGESGDVWRLVTELADWIDQDQETTPKWMPNFEALQSLS".to_vec(), b"MTSETQNSVSPPLAMPGAVAVNPRKRGRTAYVADDASSIAYTRALEERVAFLENKLAQVPTPEATTTPRETASNYSVPSGRDKNALSDVVAHVSLGNFEAPAYVGPSSGLSLALNLGEMVQATVWNKMLPDIQDGTTGNQANCINPSPRCITVEDLLAHSVKEPPSDEQGSQMLKAYTSQLHSKYPFLEPEELWKLHSERLTLAAKPTQTLTRIERFGIFKLYLVYAMGATLVQLTQRGPVLSPEALYITALQHISAARESRTVQNIEAMTLLVMFHLRSTSSHGLWYMIGLAMRTSIDLGLHRAAHEQNLDGPIVQRRRRLFWSVYSLERTIAVSLGRPLSIADNQIDVELPNTSINESPSASVIVGNDITLALVLFKLRRIESKIHHSVYRTDKTLDSLRPKLDRLHQQLKIWRNSLTDWIPTGHPDLNYALLLYNRALRLLIQPFLPILPATDPFYGLCMRAAGDICQAHKRLHQTLDYGHSFIAVQTVFVAGVTLVYGLWTQGNALWSVAVSNDIRACSLVLFVMSERAPWVRKYRDAFEVLVNAAMEKLQDSEAGLAEMASAQMRAGKAPGAADSRGVQNPDLSGNETTTRPMDSSSNQFLMSEDGGIALGEFEGAWPMVAELANWIDQDTEGGSPVWMPNFELLQSLSGTWNE".to_vec()),

        // uc30_70_80
        (b"MATFVGLSTSAGRDWTKIEKLASSMFCPLKLILMPVLLDYSLGLNDLIELTVHVGDSALLGCVFQITEEKCVTKVDWMFSSGEHAKDDYVLYYYANLSVPVGRFQNRVSLVGDILRNDGSLLLENVEEADQGTYTCEIRLEKESLVFKKAVALHVLPEEPKELTVHVGDSTQLGCVFQSTEEKRMTRVDWTFSSGEHTKEEVVLRYYPKPSVPVGYFQGWGRFQNRVTLVGDTSYNDASILLQGVKESDRGSYTCSIHLGNLTFRKTTVLRVIVKEPQTSVTPLALRPEILGGNQLVIIVGIVCGTILLLPVLILIVKRTHRNKSSALGQNRKKGSIFSGRCRGQMVKRSKAKGWEGASAGSSGGFGANSAWPPPWGRSPWSWVSLSFCCPLPAQPHLPRPGFLQHPIPWRPTLLTHLKLCGQKDGS".to_vec(), b"MFYPPKRILVPVLLSYFLGLNDLIVSSVELTVHVGDSALLGCIFQSTEEKLVTKVDWMFSSGEHFKDDYVLFYYANISVPVGRFQNRVSLVGDILHHDGSLLLQNVEEADQGNYTCEIRFKMESLVFKKAVVLHVLPEEPKELMAHVGDSTQMGCVFHSTEEKHMTRVDWMFSSGEHTKEEIVLRYYPKLKAAMGYPQNWGRFQNRVNLVGDTSHNDGSIVLHRVKESDGGSYTCSIHLGNLTVRKTTVLHVILKEPRTLVTSVTLRPEILGGNQLVIIVGVVCATILLLPVLILIVKRTYGNKSSVTSTTLVKNLENTKKANPEKHIYSSITMQEVTDEGSSGKSEATYMTMHPVWPSLRSAPTSPSDKKSDGGMPRTEQAF".to_vec())
    ];

    let cell_size = 1;
    let bg_color = Rgb([255u8, 255u8, 255u8]);
    let fg_colors = [Rgb([50u8, 50u8, 50u8]), Rgb([50u8, 50u8, 50u8]), Rgb([50u8, 50u8, 50u8])];
    let trace_color = Rgb([255u8, 0u8, 0u8]);

    for (i, img_path) in args.enumerate() {
        let q = &seqs[i].0;
        let r = &seqs[i].1;

        let r_padded = PaddedBytes::from_bytes::<AAMatrix>(r, 2048);
        let q_padded = PaddedBytes::from_bytes::<AAMatrix>(q, 2048);
        let run_gaps = Gaps { open: -11, extend: -1 };

        let mut block_aligner = Block::<true, false>::new(q.len(), r.len(), 256);
        block_aligner.align(&q_padded, &r_padded, &BLOSUM62, run_gaps, 32..=256, 0);
        let blocks = block_aligner.trace().blocks();
        let mut cigar = Cigar::new(q.len(), r.len());
        block_aligner.trace().cigar(q.len(), r.len(), &mut cigar);

        let img_width = ((r.len() + 1) * cell_size) as u32;
        let img_height = ((q.len() + 1) * cell_size) as u32;
        let fg_color = fg_colors[i];
        let mut img = RgbImage::new(img_width, img_height);

        println!("path: {}, img size: {} x {}", img_path, img_width, img_height);

        draw_filled_rect_mut(&mut img, Rect::at(0, 0).of_size(img_width, img_height), bg_color);

        for block in &blocks {
            if block.width == 0 || block.height == 0 { continue; }
            let x = (block.col * cell_size) as i32;
            let y = (block.row * cell_size) as i32;
            let width = (block.width * cell_size) as u32;
            let height = (block.height * cell_size) as u32;

            draw_filled_rect_mut(&mut img, Rect::at(x, y).of_size(width, height), fg_color);
            draw_hollow_rect_mut(&mut img, Rect::at(x, y).of_size(width, height), bg_color);
        }

        let mut x = cell_size / 2;
        let mut y = cell_size / 2;
        let vec = cigar.to_vec();

        for op_len in &vec {
            let (next_x, next_y) = match op_len.op {
                Operation::M => (x + op_len.len * cell_size, y + op_len.len * cell_size),
                Operation::I => (x, y + op_len.len * cell_size),
                _ => (x + op_len.len * cell_size, y)
            };
            draw_line_segment_mut(&mut img, (x as f32, y as f32), (next_x as f32, next_y as f32), trace_color);
            x = next_x;
            y = next_y;
        }

        let writer = BufWriter::new(File::create(img_path).unwrap());
        let encoder = PngEncoder::new_with_quality(writer, CompressionType::Best, FilterType::Sub);
        encoder.encode(img.as_raw(), img.width(), img.height(), ColorType::Rgb8).unwrap();
    }
}
