import React, { useEffect, useState } from 'react';
import Layout from '../components/Layout';
import { useRouter } from 'next/router';
import axios from 'axios';

const statusPage = () => {
  const router = useRouter();
  // Retrieve the jobHash from the URL query parameters
  const { jobHash } = router.query;

  const [isLoading, setIsLoading] = useState(false);


  useEffect(() => {
    const checkReportStatus = async () => {
      // Ensure jobHash is defined before making the API call
      if (!jobHash) {
        return;
      }

    try {
      // Make an API call to check the job status using the jobHash 
      setIsLoading(true);
      const response = await axios.get(`/api/${jobHash}`);
      const data = response.data;
      setIsLoading(false);

      // Redirect to the report
      if (data.status === 'completed') {
        window.location.href = `/reports/${jobHash}/interface/overview.html`;
      }
    } catch (error) {
      console.error('Error checking report status:', error);
      setIsLoading(false);
    }
  };
  
    // Only set up the interval if jobHash is defined
    if (jobHash) {
      const intervalId = setInterval(() => {
        checkReportStatus();
      }, 5000);
  
      return () => clearInterval(intervalId);
    }
  }, [jobHash, router]);

  return (
    <Layout>
      <h1 className="h3 mb-3 text-gray-800">  
      {isLoading ? "Fetching data..." : "The report is being generated. Please wait..."}
      </h1>
    </Layout>
  );
};

export default statusPage;