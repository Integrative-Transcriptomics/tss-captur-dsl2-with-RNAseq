import React, { useEffect, useState } from 'react';
import axios from 'axios';

export default function UploadPage() {
  // State to keep track of the motifNumber
  const [motifNumber, setMotifNumber] = useState(5);
  // State to keep track of the max file size and total file size
  const [maxFileSize, setMaxFileSize] = useState(null);
  const [maxTotalFileSize, setMaxTotalFileSize] = useState(null);

  useEffect(() => {
    // Fetch the file size limits from the API when the component mounts
    const fetchConfig = async () => {
      try {
        const response = await axios.get('/api/config');
        const config = response.data;
        setMaxFileSize(config.maxFileSize);
        setMaxTotalFileSize(config.maxTotalFileSize);
      } catch (error) {
        console.error('Error fetching config:', error);
      }
    };

    fetchConfig();
  }, []);

  // Helper function to validate file sizes
  const validateFileSizes = (fileInputs) => {
    let totalFileSize = 0;
    for (let input of fileInputs) {
      for (let file of input.files) {
        totalFileSize += file.size;
        if (file.size > maxFileSize) {
          alert(`File ${file.name} exceeds the maximum file size limit.`);
          return false;
        }
        if (totalFileSize > maxTotalFileSize) {
          alert('The total size of the files exceeds the maximum total file size limit.');
          return false;
        }
      }
    }
    return true;
  };
  const handleUpload = async (event) => {
  };

  return (
    <div>
      <h1>Upload Files for TSS-CAPTUR</h1>
      <p>Maximum size per file: {maxFileSize / 1024 / 1024} MB</p>
      <p>Maximum total file size: {maxTotalFileSize / 1024 / 1024} MB</p>
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
