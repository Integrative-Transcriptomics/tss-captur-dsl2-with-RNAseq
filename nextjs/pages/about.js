import Layout from '../components/Layout';

export default function About() {
  return (
    <Layout>
      <div>
        <h1 className="h3 mb-3 text-gray-800">About Page</h1>
        <div>
          Developed by Mathias Witte Paz
          <br />Webified by Thomas Vogel
        </div>
        <p>
          <br />To enable the characterization of transcripts in an explorative manner, we have developed TSS-Captur. TSS-Captur characterizes transcripts starting on clear TSS signals obtained from TSSpredator (Herbig, 2015) that cannot be associated to any known gene. To characterize the function of the transcript, TSS-Captur integrates tools using different approaches. The following figure describes the different steps and tools run by TSS-Captur.
          <img src="https://user-images.githubusercontent.com/29492782/119950451-2403d880-bf9b-11eb-9b3e-326408f47c53.png" style={{ width: '75%', display: 'block', margin: 'auto' }}></img>
          <br />It combines the advantages of comparative genomics and ab initio methods to account for many characteristics in the prediction. Furthermore, it uses the approach of transcriptional-signal based tools to combine TSS signals with a possible terminator. Moreover, it enables the user to explore the thermodynamic characteristics of possible identified ncRNA transcripts and visualize the secondary structure. Lastly, it facilitates the analysis of motifs within the promoter regions before each TSS signal to identify possible binding sites or other regulatory elements. All these results are combined within an interactive interface. This interface does not only provide an overview of the predicted results, but offers the possibility of getting the detailed information for each transcript and hence being able to reconstruct the predictions done by TSS-Captur.
        </p>
        <p>
          References:
          <br /><a href="https://doi.org/10.1093/nar/gkz400">CNIT <i class="fas fa-external-link-alt"></i></a>
          <br /><a href="https://doi.org/10.1186/1471-2105-2-8">QRNA <i class="fas fa-external-link-alt"></i></a>
          <br /><a href="https://doi.org/10.1016/S0022-2836(05)80360-2">BLAST <i class="fas fa-external-link-alt"></i></a>
          <br /><a href="https://doi.org/10.1186/s12859-019-2704-x">RhoTermPredict <i class="fas fa-external-link-alt"></i></a>
          <br /><a href="https://doi.org/10.1186%2Fgb-2007-8-2-r22">TransTermHP <i class="fas fa-external-link-alt"></i></a>
          <br /><a href="https://doi.org/10.1093%2Fnar%2Fgkn188">RNAfold <i class="fas fa-external-link-alt"></i></a>
          <br /><a href="https://doi.org/10.1093%2Fnar%2Fgkv416">MEME <i class="fas fa-external-link-alt"></i></a>
        </p>
      </div>
    </Layout>
  );
}