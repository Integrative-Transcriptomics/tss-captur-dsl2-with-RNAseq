import React, { useEffect, useState } from 'react';
import Layout from '../components/Layout';
import { useRouter } from 'next/router';
import axios from 'axios';

const statusPage = () => {
  // Retrieve the jobHash from the URL query parameters
  const router = useRouter();
  const { jobHash } = router.query;
  // State to keep track of the report status
  const [status, setStatus] = useState('pending');
  // State to contain an error message
  const [errorMessage, setErrorMessage] = useState('');
  // State to store the interval for status check
  const [checkInterval, setCheckInterval] = useState(5000); // Default value

  useEffect(() => {
    // Fetch the interval from the config API
    const fetchConfig = async () => {
      try {
        const response = await axios.get('/api/config');
        const config = response.data;
        setCheckInterval(config.checkInterval);
      } catch (error) {
        console.error('Error fetching config:', error);
      }
    };

    fetchConfig();
  }, []);

  useEffect(() => {
    let intervalId;

    // Check report status if the current status is 'pending'
    if (status === 'pending') {
      // Function to check the report status
      const checkReportStatus = async () => {
        // Ensure jobHash is defined before making the API call
        if (!jobHash) {
          setStatus('error'); 
          setErrorMessage('Job hash is missing.');
          return;
        }

        try {
          // Make an API call to check the job status using the jobHash 
          setStatus('fetching');
          const response = await axios.get(`/api/${jobHash}`);
          const data = response.data;

          // Handle fetched report status based on the response
          if (data.status === 'completed') {
            // Redirect to the report interface if the job is completed
            window.location.href = `/reports/${jobHash}/Interface/overview.html`;
          } else if (data.status === 'failed') {
            // Update state to 'failed' and sets the error message
            setStatus('failed');
            setErrorMessage(data.errorMessage);
          } else {
            // Keep the status as 'pending' if the job is not yet completed or failed
            setStatus('pending');
          }
        } catch (error) {
          console.error('Error checking report status:', error);
          setStatus('error'); 
          setErrorMessage('An error occurred while fetching the report status.')
        }
      };
    
      // Set an interval to periodically check the report status
      intervalId = setInterval(() => {
        checkReportStatus();
      }, checkInterval);
    }

    // Cleanup function to clear the interval on component unmount
    return () => clearInterval(intervalId);
  }, [jobHash, status, checkInterval]);

  // Render status page
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