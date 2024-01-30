import React, { useEffect, useState } from 'react';
import Layout from '../components/Layout'; // Adjust the import path as necessary
import { useRouter } from 'next/router';
import axios from 'axios';

const statusPage = () => {
  const router = useRouter();
  const { jobHash } = router.query;
  const [reportReady, setReportReady] = useState(false);

  useEffect(() => {
    const checkReportStatus = async () => {
      // Ensure jobHash is defined before making the API call
      if (!jobHash) {
        return;
      }
  
      try {
        const response = await axios.get(`/api/${jobHash}`);
        const data = response.data;

        if (data.status === 'completed') {
          setReportReady(true);
          // Redirect to the report
          window.location.href = `/reports/${jobHash}/interface/overview.html`;
        }
      } catch (error) {
        console.error('Error checking report status:', error);
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
      <h1 className="h3 mb-3 text-gray-800">The report is being generated. Please wait...</h1>
    </Layout>
  );
};

export default statusPage;