import { spawn } from 'child_process';
import path from 'path';

export default function runHandler(req, res) {
  try {
    if (req.method !== 'POST') {
      return res.status(405).end('Method Not Allowed');
    }

    const { jobHash } = req.body;

    if (!jobHash) {
      return res.status(400).end('Missing jobHash');
    }

/*     const paramsPath = path.join('../jobs', jobHash, 'params.json');
    const scriptPath = path.join('../nf', 'tss_captur.sh');

    // Start the Nextflow script asynchronously
    const script = spawn('bash', [scriptPath, '-params-file', paramsPath]);

    script.stdout.on('data', (data) => {
      console.log(`stdout: ${data}`);
    });

    script.stderr.on('data', (data) => {
      console.error(`stderr: ${data}`);
    });

    script.on('error', (error) => {
      console.error(`Error: ${error.message}`);
    });

    script.on('close', (code) => {nextjs/pages/api/upload.js
      console.log(`child process exited with code ${code}`);
    }); */

    // Return an immediate response with the status URL
    const reportUrl = `http://localhost:3000/report/${jobHash}`;
    res.status(200).send({ message: 'Analysis started successfully', reportUrl: reportUrl });
  } catch (err) {
    console.error(`Unhandled Error: ${err.message}`);
    return res.status(500).send('Unhandled Error');
  }
}
