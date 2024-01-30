import { spawn } from 'child_process';
import fs from 'fs';
import path from 'path';

export default async function runHandler(req, res) {
  try {
    if (req.method !== 'POST') {
      return res.status(405).end('Method Not Allowed');
    }

    const { jobHash } = req.body;

    if (!jobHash) {
      return res.status(400).end('Missing jobHash');
    }

    const jobDir = path.join('./public/jobs', jobHash);
    const reportDir = path.join(jobDir, 'interface');
    const paramsDir = path.join(jobDir, 'params.json');
    const scriptDir = path.join('../nf', 'tss_captur.sh');
    const placeholderDir = path.join('./static', 'placeholder.html');
    const targetDir = path.join(reportDir, 'overview.html');

    // Create placeholder page after successful run start
    try {
      await fs.promises.mkdir(reportDir, { recursive: true });
      await fs.promises.copyFile(placeholderDir, targetDir);
    } catch (err) {
      console.error('Error during report directory creation:', err);
      throw err;
    }

    /*     

    // Start the Nextflow script asynchronously
    const script = spawn('bash', [scriptDir, '-params-file', paramsDir]);

    script.stdout.on('data', (data) => {
      console.log(`stdout: ${data}`);
    });

    script.stderr.on('data', (data) => {
      console.error(`stderr: ${data}`);
    });

    script.on('error', (error) => {
      console.error(`Error: ${error.message}`);
    });

    script.on('close', (code) => {
      console.log(`child process exited with code ${code}`);
    }); */

    // Return an immediate response with the report URL
    const reportUrl = `http://localhost:3000/jobs/${jobHash}/interface/overview.html`;
    res.status(200).send({ message: 'Analysis started successfully', reportUrl: reportUrl });
  } catch (err) {
    console.error(`Unhandled Error: ${err.message}`);
    return res.status(500).send('Unhandled Error');
  }
}
