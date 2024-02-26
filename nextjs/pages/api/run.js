import { spawn } from 'child_process';
import fs from 'fs';
import path from 'path';

export default async function runHandler(req, res) {
  try {
    // Check if the request method is POST
    if (req.method !== 'POST') {
      return res.status(405).end('Method Not Allowed');
    }

    // Extract jobHash from the request body
    const { jobHash } = req.body;

    // Validate the presence of jobHash
    if (!jobHash) {
      return res.status(400).end('Missing jobHash');
    }

    // Define directories for uploads and job processing
    const uploadDir = path.resolve('../nextjs/uploads', jobHash);
    const jobDir = path.join('./public/reports', jobHash);
    const paramsDir = path.join(uploadDir, 'params.json');

    // Attempt to create the job directory
    try {
      await fs.promises.mkdir(jobDir, { recursive: true });
    } catch (err) {
      console.error('Error during job directory creation:', err);
      throw err;
    }  

    // Execute the Nextflow script asynchronously
    console.log('trying to run');
    const command = `./nextflow run ./tss_captur.nf -params-file ${paramsDir}`;
    const script = spawn('bash', ['-c', command], { cwd: '../nf' });

    // Handle script stdout stream
    script.stdout.on('data', (data) => {
      console.log(`stdout: ${data}`);
    });

    // Handle script stderr stream
    script.stderr.on('data', (data) => {
      console.error(`stderr: ${data}`);
    });

    // Handle script error event
    script.on('error', (error) => {
      console.error(`Error: ${error.message}`);
    });

    // Handle script close event
    script.on('close', (code) => {
      console.log(`child process exited with code ${code}`);
    }); 

    // Respond immediately with the report URL
    res.status(200).send({ message: 'Analysis started successfully', reportUrl: `/status?jobHash=${jobHash}` });
  } catch (err) {
    console.error(`Unhandled Error: ${err.message}`);
    res.status(500).send('Internal Server Error');  }
}
