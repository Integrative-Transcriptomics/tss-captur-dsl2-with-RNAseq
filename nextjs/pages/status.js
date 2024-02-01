import React, { useEffect, useState } from 'react';
import Layout from '../components/Layout';
import { useRouter } from 'next/router';
import axios from 'axios';

const statusPage = () => {
  const router = useRouter();
  // Retrieves the jobHash from the URL query parameters
  const { jobHash } = router.query;
  // State to keep track of the report status
  const [status, setStatus] = useState('pending');
  // State to contain an error message
  const [errorMessage, setErrorMessage] = useState('');


  useEffect(() => {
    const checkReportStatus = async () => {
      // Ensures jobHash is defined before making the API call
      if (!jobHash || (status !== 'pending' && status !== 'error')) {
        return;
      }

      try {
        // Makes an API call to check the job status using the jobHash 
        setStatus('fetching');
        const response = await axios.get(`/api/${jobHash}`);
        const data = response.data;

        // Handles fetched report status
        if (data.status === 'completed') {
          window.location.href = `/reports/${jobHash}/interface/overview.html`;
        } else if (data.status === 'failed') {
          setStatus('failed');
          setErrorMessage(data.errorMessage);
        } else {
          setStatus('pending');
        }
      } catch (error) {
        console.error('Error checking report status:', error);
        setStatus('error'); 
        setErrorMessage('An error occurred while fetching the report status.')
      }
    };
  
    const intervalId = setInterval(() => {
      checkReportStatus();
    }, 5000); // every 5 seconds

    return () => clearInterval(intervalId);
  }, [jobHash, status]);

  return (
    <Layout>
      <h1 className="h3 mb-3 text-gray-800"> 
        {status === 'fetching' && "Fetching data..."} 
        {status === 'pending' && "Job is still running. Please wait..."}
        {status === 'failed' && <>
          <div>Job failed. Error log:</div>
          <p>{errorMessage}</p>
        </>}
        {status === 'error' && errorMessage}
      </h1>
    </Layout>
  );
};

export default statusPage;