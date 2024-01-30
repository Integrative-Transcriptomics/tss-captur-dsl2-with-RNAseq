import fs from 'fs';
import path from 'path';

export default async function handler(req, res) {
  const { jobHash } = req.query;
  const reportDir = path.join('./public/reports', jobHash, 'interface', 'overview.html');

  try {
    // TODO: Catch runtime error
    if (fs.existsSync(reportDir)) {
      // If the overview.html file exists, the job is considered completed
      res.status(200).json({ status: 'completed' });
    } else {
      // If the overview.html file doesn't exist, the job is still pending or processing
      res.status(200).json({ status: 'pending' });
    }
  } catch (error) {
    console.error(`Error checking status for job ${jobHash}:`, error);
    res.status(500).json({ error: 'Internal server error' });
  }
}