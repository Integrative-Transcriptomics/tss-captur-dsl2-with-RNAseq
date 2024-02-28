import fs from 'fs';
import path from 'path';

// Handler function to check the status of a job based on its hash
export default async function handler(req, res) {
  // Extract jobHash from the request query
  const { jobHash } = req.query;
  // Define directories
  const paramsDir = path.join('./uploads', jobHash, 'params.json');
  const reportDir = path.join('./public/reports', jobHash);
  const overviewDir = path.join(reportDir, 'Interface', 'overview.html');
  const errorDir = path.join(reportDir, 'error.log');

  try {
    // Check if the job's params.json file exists
    if (fs.existsSync(paramsDir)) {
      // Check if the overview.html exists to determine if the job is completed
      if (fs.existsSync(overviewDir)) {
        // If the overview.html file exists, the job is considered completed
        res.status(200).json({ status: 'completed' });
      } else if (fs.existsSync(errorDir)) {
        // If the error.log file exists, read its content and return it with the failed status
        const err = fs.readFileSync(errorDir, 'utf8');
        res.status(200).json({ status: 'failed', errorMessage: err });
      } else {
        // If none of the files exist, the job is still pending
        res.status(200).json({ status: 'pending' });
      }
    } else {
      // If the job does not exist, return with the failed status
      const err = 'Job does not exist';
      res.status(200).json({ status: 'error', errorMessage: err });
    }
  } catch (error) {
    console.error(`Error checking status for job ${jobHash}:`, error);
    res.status(500).send('Internal Server Error');
    }
}