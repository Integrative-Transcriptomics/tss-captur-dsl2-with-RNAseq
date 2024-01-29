import React, { useEffect, useState } from 'react';

export default function UploadPage() {
  // State to keep track of the motifNumber
  const [motifNumber, setMotifNumber] = useState(5);
  const handleUpload = async (event) => {
  };

  return (
    <div>
      <form onSubmit={handleUpload} encType="multipart/form-data">
        <div>
          <label htmlFor="masterTable">MasterTable.tsv:</label>
          <input type="file" id="masterTable" name="masterTable" accept=".tsv" required />
        </div>
        <div>
          <label htmlFor="genomeFolder">Genome Files:</label>
          <input type="file" id="genomeFolder" name="genome" accept=".fa,.fna,.fasta,.frn,.faa,.ffn" multiple required />
        </div>
        <div>
          <label htmlFor="gffFolder">GFF Files:</label>
          <input type="file" id="gffFolder" name="gff" accept=".gff" multiple required />
        </div>
        <div>
          <label htmlFor="motifNumber">Motif Number: {motifNumber}</label>
          <br />
          <input type="range" id="motifNumber" name="motifNumber" min="1" max="99" value={motifNumber} onChange={e => setMotifNumber(e.target.value)} required />
        </div>
      </form>
    </div>
  );
}
