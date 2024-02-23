import React, { useEffect, useState } from 'react';
import axios from 'axios';
import crypto from 'crypto';
import Layout from '../components/Layout';

export default function UploadPage() {
  // State to keep track of the motifNumber
  const [motifNumber, setMotifNumber] = useState(5);
  // State to keep track of the max file size and total file size
  const [maxFileSize, setMaxFileSize] = useState(null);
  const [maxTotalFileSize, setMaxTotalFileSize] = useState(null);
  // State to keep track of the upload status
  const [startedUpload, setStartedUpload] = useState(false);
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
      console.log(jobHash);
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
    formData.append('motifNumber', motifNumber);
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
      if (!startedUpload) {
        setStartedUpload(true);
        const response = await axios.post('/api/upload', formData, options);
        if (response.status === 200) {
          console.log('Files uploaded successfully');
          const data = response.data;
          setJobHash(data.jobHash);
        } else {
          console.error('Error uploading files');
          setStartedUpload(false);
        }
      }
    } catch (error) {
      console.error('Error uploading files:', error);
      setStartedUpload(false);
    }
  };

  // Helper function to start upload
  const handleUpload = async (event) => {
    try {
      event.preventDefault();
      const formData = new FormData();
      const fileInputs = document.querySelectorAll('input[type="file"]');

      if (validateFileSizes(fileInputs)) {
        await appendFilesAndHashes(fileInputs, formData);
        await uploadFiles(formData);
      }
    } catch (error) {
      console.error('Error during the upload process:', error);
      alert('An unexpected error occurred during the upload process.');
      setStartedUpload(false);
    }
  };

  // Helper function to start the analysis after a successful upload
  const handleAnalysisStart = async () => {
    try {
      const runResponse = await axios.post('/api/run', { jobHash: jobHash });
      if (runResponse.status === 200) {
        const data = runResponse.data;
        const baseUrl = window.location.origin;
        const fullReportUrl = baseUrl + data.reportUrl;
        setReportUrl(fullReportUrl);

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
    <Layout>
      <div>
        <h1 className="h3 mb-3 text-gray-800">Upload Files for TSS-CAPTUR</h1>
        <div className="d-flex flex-column mb-4" style={{ width: '300px' }}>
          <p>Maximum size per file: {maxFileSize ? (maxFileSize / 1024 / 1024).toFixed(2) : 'Loading...'} MB</p>
          <p>Maximum total file size: {maxTotalFileSize ? (maxTotalFileSize / 1024 / 1024).toFixed(2) : 'Loading...'} MB</p>
        </div>
        <form onSubmit={handleUpload} encType="multipart/form-data">
          <div className="form-group">
            <div className="d-flex flex-column">
              <label htmlFor="masterTable">MasterTable.tsv:</label>
              <input type="file" id="masterTable" name="masterTable" accept=".tsv" required />
            </div>
          </div>
          <div className="form-group">
            <div className="d-flex flex-column">
              <label htmlFor="genomeFolder">Genome Files:</label>
              <input type="file" id="genomeFolder" name="genome" accept=".fa,.fna,.fasta,.frn,.faa,.ffn" multiple required />
            </div>
          </div>
          <div className="form-group">
            <div className="d-flex flex-column">
              <label htmlFor="gffFolder">GFF Files:</label>
              <input type="file" id="gffFolder" name="gff" accept=".gff" multiple required />
            </div>
          </div>
          <div className="form-group d-flex flex-column" style={{ width: '300px' }}>
            <label htmlFor="motifNumber" className="mb-0 mr-2">Max motif numbers: {motifNumber}</label>
            <input type="range" id="motifNumber" name="motifNumber" min="1" max="20" value={motifNumber} onChange={e => setMotifNumber(e.target.value)} required />
          </div>
          <button type="submit" className="btn btn-primary" disabled={startedUpload}>Start Upload</button>
        </form>
        {startedUpload && <p className="mb-4">Upload progress: {uploadProgress}%</p>}
        {reportUrl && <div>Check the status of the report at: <a href={reportUrl}>{reportUrl}</a></div>}
      </div>
    </Layout>
  );
}