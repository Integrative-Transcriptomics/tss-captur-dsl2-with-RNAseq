import React, { useEffect, useState } from 'react';
import axios from 'axios';
import crypto from 'crypto';

export default function UploadPage() {
  // State to keep track of the motifNumber
  const [motifNumber, setMotifNumber] = useState(5);
  // State to keep track of the max file size and total file size
  const [maxFileSize, setMaxFileSize] = useState(null);
  const [maxTotalFileSize, setMaxTotalFileSize] = useState(null);
  // State to keep track of the upload status
  const [isUploading, setIsUploading] = useState(false);
  // State to keep track of the upload progress
  const [uploadProgress, setUploadProgress] = useState(0);
  // State to keep track of the jobHash
  const [jobHash, setJobHash] = useState(null);
  // State to keep track of the report URL
  const [reportUrl, setReportUrl] = useState(null);

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

  useEffect(() => {
    // Start the analysis when jobHash changes and is not null
    if (jobHash) {
      handleAnalysisStart();
    }
  }, [jobHash]);

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

  // Helper function to append files to FormData and generate SHA-256 hash for each file
  const appendFilesAndHashes = async (fileInputs, formData) => {
    const fileHashes = {};

    for (let input of fileInputs) {
      for (let file of input.files) {
        formData.append(input.name, file);
        const hash = crypto.createHash('sha256');
        const buffer = await file.arrayBuffer();
        hash.update(new Uint8Array(buffer));
        const hashFile = hash.digest('hex');
        fileHashes[file.name] = hashFile;
      };
    };
    formData.append('fileHashes', JSON.stringify(fileHashes));
  };

  // Helper function to handle the upload process
  const uploadFiles = async (formData) => {
    const options = {
      headers: { "Content-Type": "multipart/form-data" },
      onUploadProgress: (progressEvent) => {
        const percentage = (progressEvent.loaded * 100) / progressEvent.total;
        setUploadProgress(+percentage.toFixed(2));
      },
    };

    try {
      const response = await axios.post('/api/upload', formData, options);
      if (response.status === 200) {
        console.log('Files uploaded successfully');
        const data = response.data;
        setIsUploading(false);
        setJobHash(data.jobHash);
      } else {
        console.error('Error uploading files');
        setIsUploading(false);
      }
    } catch (error) {
      console.error('Error uploading files:', error);
      setIsUploading(false);
    }
  };

  // Helper function to start upload
  const handleUpload = async (event) => {
    try {
      event.preventDefault();
      setIsUploading(true);
      const formData = new FormData();
      const fileInputs = document.querySelectorAll('input[type="file"]');

      if (!validateFileSizes(fileInputs)) {
        setIsUploading(false);
        return;
      }

      await appendFilesAndHashes(fileInputs, formData);
      await uploadFiles(formData);
    } catch (error) {
      console.error('Error during the upload process:', error);
      alert('An unexpected error occurred during the upload process.');
      setIsUploading(false);
    }
  };

  // Helper function to start the analysis after a successful upload
  const handleAnalysisStart = async () => {
    try {
      const runResponse = await axios.post('/api/run', { jobHash: jobHash });
      if (runResponse.status === 200) {
        const data = runResponse.data;
        setReportUrl(data.reportUrl);
        console.log('Run API called successfully');
      } else {
        console.error('Error calling Run API');
        alert('Analysis could not be started. Please try again.');
      }
    } catch (error) {
      console.error('Error calling Run API:', error);
      alert('An error occurred while starting the analysis.');
    }
  };

  return (
    <div>
      <h1>Upload Files for TSS-CAPTUR</h1>
      <p>Maximum size per file: {maxFileSize / 1024 / 1024} MB</p>
      <p>Maximum total file size: {maxTotalFileSize / 1024 / 1024} MB</p>
      <p>Upload progress: {uploadProgress}%</p>
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
        <button type="submit" disabled={isUploading}>Start</button>
      </form>
      {reportUrl && <div>Report will be available at: {reportUrl}</div>}
    </div>
  );
}
