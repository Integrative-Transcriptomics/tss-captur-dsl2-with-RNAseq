import fs from 'fs';
import path from 'path';

export default async function handler(req, res) {
  const { jobHash } = req.query;
  const reportDir = path.join('./public/reports', jobHash);
  const overviewDir = path.join(reportDir, 'interface', 'overview.html')
  const errorDir = path.join(reportDir, 'error.log');

  try {
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
  } catch (error) {
    console.error(`Error checking status for job ${jobHash}:`, error);
    res.status(500).json({ error: 'Internal server error' });
  }
}