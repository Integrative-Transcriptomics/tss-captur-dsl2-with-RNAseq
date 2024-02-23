import Layout from '../components/Layout';

export default function Example() {
  return (
    <Layout>
      <div>
        <h1 className="h3 mb-3 text-gray-800">Example Data Page</h1>
        <p>
          Campylobacter jejuni, a Gram-negative bacterium, is known to be one of the leading causes of bacterial gastroenteritis in humans. Its adaptability to a wide range of hosts and conditions, as well as its comparatively small genome of about 1.6 Mb, make it a good choice for the study of prokaryotic transcriptional mechanisms.
          The following dataset includes four C. jejuni strains: <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=Info&id=192222">NCTC11168 <i class="fas fa-external-link-alt"></i></a>, <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?mode=info&id=354242">81-176 <i class="fas fa-external-link-alt"></i></a>, <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=407148">81116 <i class="fas fa-external-link-alt"></i></a> and <a href="https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=195099">RM1221 <i class="fas fa-external-link-alt"></i></a>. These strains were subjected to a high resolution transcriptome analysis by <a href="https://doi.org/10.1371/journal.pgen.1003495">Dugar et al. (2013) <i class="fas fa-external-link-alt"></i></a>.
          <br />The dataset is available for download here: <a href="/example/input/dataset.zip"><i class="fas fa-download"></i></a>
          <br />The corresponding report from TSS-Captur is available here: <a href="/example/report/Interface/overview.html">[Click here]</a>
        </p>
      </div>
    </Layout>
  );
}