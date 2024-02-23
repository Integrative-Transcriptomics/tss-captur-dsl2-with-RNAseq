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

    const uploadDir = path.resolve('../nextjs/uploads', jobHash);
    const jobDir = path.join('./public/reports', jobHash);
    const paramsDir = path.join(uploadDir, 'params.json');

    // Create jobDir after successful upload
    try {
      await fs.promises.mkdir(jobDir, { recursive: true });
    } catch (err) {
      console.error('Error during job directory creation:', err);
      throw err;
    }  

    // Start the Nextflow script asynchronously
    console.log('trying to run');
    const command = `./nextflow run ./tss_captur.nf -params-file ${paramsDir}`;
    const script = spawn('bash', ['-c', command], { cwd: '../nf' });

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
    }); 

    // Return an immediate response with the report URL
    res.status(200).send({ message: 'Analysis started successfully', reportUrl: `/status?jobHash=${jobHash}` });
  } catch (err) {
    console.error(`Unhandled Error: ${err.message}`);
    return res.status(500).send('Unhandled Error');
  }
}
