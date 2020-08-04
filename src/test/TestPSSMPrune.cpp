//  main.cpp
//  forautocompl
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include <iostream>
#include "Parameters.h"
#include "StripedSmithWaterman.h"
#include "MsaFilter.h"
#include "PSSMCalculator.h"
#include "Sequence.h"
#include "SubstitutionMatrix.h"
#include "MultipleAlignment.h"
#include <string.h>

const char* binary_name = "test_pssmprune";

int main (int, const char**) {
    Parameters& par = Parameters::getInstance();
    SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, -0.0);
    std::cout << "Subustitution matrix:";
    SubstitutionMatrix::print(subMat.subMatrix,subMat.num2aa,subMat.alphabetSize);
    //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
    const char *seqs[1001];
    int counter = 0;


//    seqs[counter++] = "DRETELSEEEPGMEKKMEELGGEDELVDPLTTVREQCEQLEKCVKARERLELCDQRVSSRSQTEEDCTEELFDFLHARDHCVAHKLFNSLK";
//    seqs[counter++] = "-----------------EELGEEEDLVDPLTTVREQCEQLEKCVKARERLELCDQRVSSRSQTEEDCTEELFDFLHARDHCVAHKLFNSLK";
//    seqs[counter++] = "-------------EEEEDELQSEEELVDPLTTVREQCEQLEKCVKARERLELCDERVSSRSHTEEDCTEELFDFLHARDHCVAHKLFNNLK";
//    seqs[counter++] = "------------------------ELVDPLTTIREHCEQTEKCVKARERLELCDARVSSRSHTEEQCTEELFDFLHARDHCVAHKLFNKLK";
//    seqs[counter++] = "----------------------EEELVDPLTTVREQCEQLEKCVKARERLELCDKRVSSR-------------------------------";

//    seqs[counter++] = "FKLVTTHKRTWYDKTMTNVASLPKLGVNPGLAFPTKELFLRYLDLAEVPKEARKGETMQQLGVMMEAVKVKNDKALVAVTSVGGSSIHWTPTGPDLTLSWITDKKLEEEDLPIPLSLALGTVGMPGLTAYFGLLEICGVKGGETVMVNAAAGAVGSVVGQIAKLKGCKVVGAVGSDEKVAYLQKLGFDVVFNYKTVESLEETLKKASPDGYDCYFDNVGGEFSNTVIGQMKKFGRIAICGAISTYNRTGPLPPGPPPEIVIYQELRMEAFVVYRWQGDARQKALKDLLKWVLEGEKNGMQIYMYFEAFPIIEKMLTGGKALKDNKVGAI";
//    seqs[counter++] = "---------------------------------------------------------------------------------------------------------------------AMSVLGMTGLTAYFGMTEVGQPKPGDTVVVSGAAGATGMVAGQIAKIKGAKVVGLAGSAEKCAFLRELGFDAAINYKD-KDWKKQLKDATPEYIDVFFDNTGGEILDACLARAARDARFAICGAISQYNSA--KPQGPASFMVISQRVTMKGFIVFDY-AKKYPIALKDLSEWLTQG---------------LVHEILQGGKVLEDPK----";
//    seqs[counter++] = "----------------------------------------------------------------------------VAVVPTQIATLSWSPKHGGVASSITVGKELLDDTWPIPLSLALGTVGMPGLTAYFGLLEICGVKGGETVMVNAAAGAVGSVVGQIAKLKGCKVVGAVGSDEKVAYLQKLGFDVVFNYKTVESLEETLKKASPDGYDCYFDNVGGEFSNTVIGQMKKFGRIAICGAISTYNRTGPLPPGPPPEIVIYQELRMEAFVVYRWQGDARQKALKDLLKWVLEGEEKYIKNQEY--AFP---KFAMGG-----------";
//    seqs[counter++] = "-----------------------------GLAFPTPLVIMEICEW-DLPFKKKKGEDLDPVNLYMSAERVKVCYNATSVPVIAYGRRVWKIKGT------LEHRDIEKESVALPLSLMMGAMGMPGMAAYFGFLELCQPKAGETVVVSGAAGAVGSLVGQIAKIKGCKVIGFAGTEEKCKWIEELGFDFAYNYKKTD-VDTALKEAAPDSVDCYFDNVGGMFTVKVLTHMKTFGRVSICGSISNYNDTS-VPTGPLPFFIMKSQLRVEGFLVLRWYS-RWEEGETAMLQWIKEGEGHGM--------IPEKEKFLKG------------";
//    seqs[counter++] = "------------------------------------------------------------------------------------SMVGWT-HGPITT----PGDELEYDRILVPLEHFMGLLGLPGATAYHGFMELIDPKAGETVLISGAAGSVGSLVGQIAKAEGLKVIGVTGSDEKKDWLNELGFDAVINYKT-ENLEEMVAKYAPDGIDRHFENVGGKVLEVALDNMKIGGKITLCGLISTYNDTEP-------------------------------------------------------------------------------";

    seqs[counter++] = "EEQLRVKMKQMCPELSVTYPENCMLRIEMYSPMVGVEIKYLWGDPAQESEKEGENPGTSSSEVAKWARDKEKEVRNEKCQHLVLASAVDLMGEEAKESKADSISHHCPFTNKCYASKHQRLQQYLHVHPMSKGSAHVCCKKFSRVKDHRRIQSQHHVRQEQPETFTFKGESNCCTYLRTHKAGVQRTPYCCRETTLHTASDQRRFNHKRQTRGPVKSCHISKYSRNCLVDSLEKTRQYHTETGFPVKHLFLEKSSAPLRCHKACEMTFCTHRRSEHKPTCSSSFTNRDHRARFVRLNKHNSGLEPDMVHLTVNSSSSQSSCSYVCDTARKWLKKSPTTRFLKEQYNHSCSRPLKFQGKHSHNDHNICRFRPTGLKQVFYTCECHSQQDLGEIKFLTKHTPHVGEKRDCATSRSHKCFKLKAQETHKPHGATYFEKKFGCRHEQKISLCRDSLNRKTHTIRCGYQFDQHNRLTIVLSKEHKQDPKTCNPHDPAPAHNKKVLSHLCDHCFLMSTHSTSFHSSTSKGWFAKGEFQLCQKKSTKKRDHEVKKLFVCCGELKIKTNHLDSFSWHKPCISKEHQTLPGKCYHTFRFKSLCITSAEKRHRHPETPDHCHGKGFKQKITVKCSTSRTQAKRHLHYEPPTSKCQNIGDTQSGSCKKQTEPHLCGFHIHVCEYKKEFTNHRGVLSRHRTTESYRKEQPSTPGKIKRSQCKVHLVVCLSNNAEPFGDLHAIALREFEPNEIRTPHTTPRYCSEDKHRGPLALCFKTTPHKKRLVEHRCKFDFQTDCHQNKTLRETYERHPVCCSGYPKENHQLPHSVQSRGHLENIRGDFCQPTCDKYIKTRTSGAKQTHTLRPGEMKHSQKCYQISNDLCQKTFECTIRQPECFHDVKKTCFGEKHRGIYKSSRIHSPKQHDRHICCSLPSSNFLVKLHSLGKPQYGPKSVNSSAHLQVPTQELPAFHEGPVEHVCADYDIPGTSYSDGIFNSTCGEKCAEARCFIAHRDIHRKVESKLGRCTYFQCHRVTKKISQRSHHDEGTKVTRPLPSYYKCGLCERSFVEKSALSRHQRVHKNGQSSATESQQEPWDDVNEPLMYVKLVNKEHQIEPELYSSLISPVPHNPAPLSLIQSRIFTKVPYPEGCCKSSTCGRSEYWHVHMFGFTEKQYPCKRTSSKGFHQLKESAVITHHPSYSGLGSCFKPERKTCQVKISCLHEKRQDIGTACTHKTVGICVRFHSTKKHYKGPPTCECRDSYTKEKNSVVFHRGDTMKKHRVELSCGTYSCFKRQSTGDPHQRMHTCPAVDNNSYRGRCSNQDLEASFYSEGWEKLDIYGLTYQMLAKPESPEELIVIGAWIQCA";
    seqs[counter++] = "------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------KSFSKHASC------SNHQRVHTKRNSGYLKEFSCKFVSRHCQSYVAEHTG------------KKH--ECGECGKSFSQKR---SLIQHQRFHTGEKPCLSEKHRSFYESHKVHHQRVHERPYQCGECGKSFSQKVHTGARPECFGKSCSHGSHLHENIRK-----EG-YEHLGTFKDCKGTGRCYGSSHLICGE-CGKSFSSNGHLRSHQRVHT--GERPY-ECGECGKSFSHKRHHQRIHTGERP----YKCGECGKSFIEKGHLRNHQRVH---------TAERPFNQCG----CGKLEHKGYKVG--------SPV----------------RFLYREFCHKGCLGKHGLFKHAKRIHNGEKPYACEACQKFFRHKYQ---LIAHQRVHTG------------ERPYEC-----NDCGKSFTHSTFHVHKRIHTGEKPYE---CSECGKSFAESSSLTKHRRVHTGEKPFFKKNC--YSCEKHKKLGLVHQKVHT-----------------------------------------------------------";
    seqs[counter++] = "EEMSRMKMN---PQ------KTCCLQYELYIPMESEEVSLKW--EARPGNSEGEGTMPKGEEDVSVANDRAKEQMSHKCE--VLQSPVKKVEEKSDANIHDSLGYF--HSGECSCKSHVRMQAVKQYHAQSHHCPH-CKKSFVQRSDFIKHQRTHTGLQYREEAFK-------CKVPERHQ---QRTSTCVHNTRPYTCLDQKTFNHRRTHTGPYVSCRCSKFIQN---SDLVKHLRTHTE-GKPYE--------CPLCVKRFAESSLMKHKRTHHRTSCSRSFTHNSDLTAHMR--KHTSVGTVSNDLLVFPSRNPENSDLYSCVCTQKWRKKASSRQKSHKQSLFECY-------PNTHNIHKLCNFTHTG-ERPYQCAECHQKSDL--VKHLRTHT---GEKSHCDKKFTERSALAKHQRTHKPYKCSDCGKEF----TQRSNLI---LHQRIHTGERPYKF--HQRDTCKSLTHLQVD--CVNKSDPHPIHSKPTCTLACFSDKTMFLSHNSHHESKKKLWNG-SKFQCCKKGFTQKSD-LVKHIRVHTGEKPFKCLLCKSFSQNSD-LHKHWRIHTGECCSYFDTKEPTFTHLGKRHKITERPHKCSQKGFIQKSAL-TKHSRT----HTGEKPPCTQCGKSFIQNSDLVKHQREKPTGTFYSHCSEHCKNL--ERRKVRTHSGEKPYRCPQCEKTF-IQSSDLVKHLV---EVNPGPHAEETNIANRARPITEDHFLRSLPYPCTECGHQRPALLHLRTHTKEKRYPCNECKSFFQTSDLVKHLRTHTGERPYLKNCGPVSFDHEKNIHQCQRTHT-GERPYTCS-QCDKGFIQRSALTKRTHTEKPYKCEQCQKCFIQNSDLVKHQRIHTEKGPYHCPDCDKR-FTE-GSSLIKHQRIHRIKSYPCGVCSFSQSSNLLLKCHNESLQQASESTGELPFDVPVPY----PHHDDTVAYPIVSSEGEGARIAATSFKNDCGKCFAHRSVLIKHVRIH------TGERPYSQF--ILDSVSRRCKYCEGTK-TRTQPHHYKCGLCERSFVEKSALSRHQRVHKNWQESYQLMVQNPKNDDEESVSPVSLI---PHLEPHKVANVY---------LQISESQIFPKIELLEPTSSPPYSGRCSYLHHGCWRFKEQPYTCKE-CGKSFSQ--SSALVKHVRIHTGCSTCGKAPQTTRKIKIKCFCHEDGHRSSVVKHS------RTHTGERPYK----CECTKGFV-QKSDLVKH----MRTHTGEKPYGCNCCDRHVFRHASCQRSTTGPQMDNSSYYNSARLFYTEKEDYEQ--DACQPWGIPSGLL--------LEIIGSW----";
    seqs[counter++] = "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------CHQ-KHKTDSAQKH-ICC-----EKH-FRHSSALRVHQRIHTGEYQCNECGKYAQSNFNHHQRTHTEKPYECKECGKSFCVKSNLTEKPFKCTFFQRSQFADHQRTHTGEKEHGCYKKGEQEKNKHYEKPHKCTECGKSFCYKIHQVGKKPYDCNQLKVKSHLLGETNCKCKIHTGEKPYECSECGKSFSMKSDLVVHTHTGEKCNECECNECGKSFHMKSTLVIHQRTHFKCNEKNSLVSEQCHHLTGKKP-------HVCNECGKAFSMKSNLTDHQRTH---------SKEKPY----ECFECQKTVHKTHTLSRHRTHTGEKPYKCNECGKSFYHQRIHTEKPY--GC---KECGKAFFQKSHLTKHQEKPYECKECKKTFF---QKSHLIEHQRTHTG---THEPGKRCSKTALSCIFKHQ-------------YCERTHTEEKPYK----CECDKSFCMKSHLTVHQRTHTGKNPFGCSECG--KTFYVKSNLISHQRTHT-----------------------------------------------------------";
    seqs[counter++] = "---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ERTPLYCKGNEKGTHLNISFLKTIRHRRIHTGEYTCPECGKSF-SRSSNLVTHQRIHTGEKPYTCEFWSSSRHRCHPITEC--KRKYCFHGRSYKCVSERHEHLLGRPKIYRGREATQ-TCCGIYSSQLVTRIHTGERPYSCPKGFTSSSHLVTHRR----IHEGTRPYSCPE----------------CGKGFTSSSQLVSHGRIH------TGKKPY-TCQECGKNFRWSSHLIIHRKITGERPYCTEC--CGKSFIQSSLLNKHQRIHESGGPYRCTECGKDFIYSSQLVTHQRIHTGERPYTPECGKSFTRSFSLN------IHHRIHTEKPYTCHECGKSFIQRSELNKHQRVHTGERPYRCPE-CGKDF--IYSSQLVTHRRIHTG------PERYSCPVCGKSFSGSSQMI--------------THQSHRVGEAPYTCECGKSFSWHSNLITHQRIHTGEKHNQSTSRNTTGCFKHLSCPTPREKEHSCLIVETSRYEGKGSS--------------------------------------------";
    seqs[counter++] = "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------TCNKQVCFNTHGLKRLQTLREGPVSHGCSCTQISRKLCEKTFEC----GEC----GKT-FWEK-SNLTQHQRTHTEKPYECTECCQKPHLTNHQRTHTGEKP-YECKQCGKTFCVKSNLTEHQRTHTGEKPYEC----------------NACGKSFCHRSALTVHQRTHKKGSNKSGNCVFFCC-----IVHQRTHTGE-----KP----YKCNECGKTFCEKSALTKHQRTHTGNACGKTFSQR----SVLTKHQRIHRLSKINEICPKERS-----VSCNRGIINITSETHTGERPY--------ECDECEKTFFHKSSLTVPQRLCECSSHQQFHELRNT---TTAQKCKGKGKCSSPVHKEHTLKCDYSHGQEESMKESCRVESTSLCTREYTQEKIYE---CSECRKMFS-GKSDLLNHQRTHTGEKSYECVCSQTEGFHKKTKTLRH-KLHIKPFECNKCGKAFCQKSQL-----------------------------------------";
    seqs[counter++] = "--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------KCEFRHRGTHGHGLHQMIKSIY-ECGDCGKSFSYKSNLIEHQRVTRERP--YECGECRKSFRQSSSLFRHQRVHTGEKPYECESCGKTFRQIFNLIRHERVMSQECKLPHRIHSG-ERPYECNECGKSLIRHRVHTERPY-EVKCSCSGCFRKQHQRVHT---GERPYECGE-CGKSF--TRKSDLIQHQRIHTGTRECGRESSFSCQKHVRCIEERKGQSASLIQHQ------RVHTGEKPYE---CSECGKSFSQSSSLIQHQRGHT--GERPEYGKCNEPCFTHKSDLIQHQRVHT-----GREPYGRVTFEDVTVYFSSEEWDLLD----------------------------";
    char ** seqsCpy = new char*[counter];
    for (int k = 0; k < counter; ++k) {
        seqsCpy[k] = MultipleAlignment::initX(strlen(seqs[0]));
        for (size_t pos = 0; pos < strlen(seqs[0]); ++pos) {
//            seqs[k][pos] = (seqs[k][pos] == '-') ? MultipleAlignment::GAP : subMat.aa2num[(int) seqs[k][pos]];
            seqsCpy[k][pos] = (seqs[k][pos] == '-') ? MultipleAlignment::GAP : static_cast<int>(subMat.aa2num[(int) seqs[k][pos]]);
        }
    }

    MultipleAlignment::MSAResult res(strlen(seqs[0]), strlen(seqs[0]), counter, seqsCpy);
    MultipleAlignment::print(res, &subMat);

    MsaFilter msaFilter(10000, counter, &subMat, par.gapOpen.aminoacids, par.gapExtend.aminoacids);
    msaFilter.pruneAlignment((char**)res.msaSequence, res.setSize, res.centerLength);

    std::cout <<"Pruned MSA" << std::endl;
    for(int k = 0; k < (int)res.setSize; k++){
        //printf("k=%.3d ", k);
        for(size_t pos = 0; pos < res.centerLength; pos++){
            char aa = res.msaSequence[k][pos];
            printf("%c", (aa < MultipleAlignment::NAA) ? subMat.num2aa[(int)aa] : '-' );
        }
        printf("\n");
    }
    std::vector<Matcher::result_t> empty;
    size_t filterSetSize = msaFilter.filter(res, empty, 0, 0, -20.0f, 90, 100);
    std::cout << "Filtered:" << filterSetSize << std::endl;
//    for(size_t k = 0; k < res.setSize; k++){
//        std::cout << "k=" << k << "\t" << (int)filterResult.keep[k] << std::endl;
//    }
    std::cout <<"Filtered MSA" << std::endl;
    for(size_t k = 0; k < filterSetSize; k++){
        printf("k=%.3zu ", k);
        for (size_t pos = 0; pos < res.centerLength; pos++) {
            char aa = res.msaSequence[k][pos];
            printf("%c", (aa < MultipleAlignment::NAA) ? subMat.num2aa[(int) aa] : '-');
        }
        printf("\n");
    }


    //seqSet.push_back(s5);
//    PSSMCalculator pssm(&subMat, counter, 1.0, 1.5);
//    pssm.computePSSMFromMSA(filterResult.setSize, res.centerLength, filterResult.filteredMsaSequence, false);
//    //pssm.printProfile(res.centerLength);
//    pssm.printPSSM(res.centerLength);
    for (int k = 0; k < counter; ++k) {
        free(seqsCpy[k]);
    }
    delete [] seqsCpy;
    return 0;
}

//PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
//                     ALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL
