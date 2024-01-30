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

    const uploadDir = path.join('./uploads', jobHash);
    const jobDir = path.join('./public/reports', jobHash);
    const paramsDir = path.join(uploadDir, jobHash, 'params.json');
    const scriptDir = path.join('../nf', 'tss_captur.sh');

    // Create jobDir after successful upload
    try {
      await fs.promises.mkdir(jobDir, { recursive: true });
    } catch (err) {
      console.error('Error during job directory creation:', err);
      throw err;
    }

    // TODO: Remove debug report
    try {
      const interfaceDir = path.join(jobDir, "interface");
      await fs.promises.mkdir(interfaceDir, { recursive: true });
      const overviewDir = path.join(interfaceDir, "overview.html");
      await fs.promises.copyFile("./static/overview.html", overviewDir);
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
    res.status(200).send({ message: 'Analysis started successfully', reportUrl: `/status?jobHash=${jobHash}` });
  } catch (err) {
    console.error(`Unhandled Error: ${err.message}`);
    return res.status(500).send('Unhandled Error');
  }
}
