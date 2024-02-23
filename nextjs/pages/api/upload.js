import { IncomingForm } from 'formidable';
import fs from 'fs';
import path from 'path';
import crypto from 'crypto';
import { MAX_FILE_SIZE, MAX_TOTAL_FILE_SIZE, calculateUploadTimeout } from './config';

export const config = {
  api: {
    bodyParser: false,
  },
};

const UPLOADS_DIR = path.resolve('./uploads');
const REPORT_DIR = path.resolve('./public/reports')

export default async function uploadHandler(req, res) {
  let jobHash;
  let jobUploadDir;
  
  try {
    if (req.method !== 'POST') {
      return res.status(405).end('Method Not Allowed');
    }

    // Generate a unique hash for the job
    do {
      jobHash = crypto.randomBytes(8).toString('hex');
      jobUploadDir = path.join(UPLOADS_DIR, jobHash);
    } while (fs.existsSync(jobUploadDir));

    console.log(jobHash);
    // Create job specific directories
    const genomeDir = path.join(jobUploadDir, 'genome');
    const gffDir = path.join(jobUploadDir, 'gff');

    await fs.promises.mkdir(jobUploadDir, { recursive: true });
    await fs.promises.mkdir(genomeDir, { recursive: true });
    await fs.promises.mkdir(gffDir, { recursive: true });

    // Form options
    const form = new IncomingForm({ hashAlgorithm: 'sha256' });
    form.uploadDir = jobUploadDir;
    form.keepExtensions = true;
    form.multiples = true;
    form.maxFileSize = MAX_FILE_SIZE;
    form.maxTotalFileSize = MAX_TOTAL_FILE_SIZE;

    let masterTableFile;
    let motifNumber;
    let clientHashes;

    const uploadTimeout = calculateUploadTimeout();

    await new Promise((resolve, reject) => {
      const timer = setTimeout(() => reject(new Error('Upload timeout exceeded')), uploadTimeout);

      form.on('error', reject);

      // Rename & move files to appropriate destinations
      form.on('fileBegin', async function(field, file) {
        let newPath;
        if (field === 'masterTable') {
          newPath = path.join(jobUploadDir, file.originalFilename);
          masterTableFile = file.originalFilename;
        } else {
          newPath = path.join(jobUploadDir, field, file.originalFilename);
        }
        file.filepath = newPath;
      });

      // Resolves after getting motifNumber and comparing hashes
      form.parse(req, (err, fields, files) => {
        clearTimeout(timer);
        if (err) return reject(err);

        motifNumber = parseInt(fields.motifNumber, 10);
        clientHashes = JSON.parse(fields.fileHashes);
        // Compare the hashes
        for (let field in files) {
          for (let file of files[field]) {
            const formHash = file.hash;
            const clientHash = clientHashes[file.originalFilename];
            if (formHash !== clientHash) return reject(new Error('Hash mismatch detected'));
          }
        }
        resolve();
      });
    });
    console.log('Creating params');
    // Create params.json after successful upload
    const jobReportDir = path.join(REPORT_DIR, jobHash);
    const tableDir = path.join(jobUploadDir, masterTableFile);
    const params = {
      inputTable: tableDir,
      inputGenomes: genomeDir,
      inputGFFs: gffDir,
      motifNumber: motifNumber,
      outputDir: jobReportDir,
    };
    await fs.promises.writeFile(path.join(jobUploadDir, 'params.json'), JSON.stringify(params));

    // Respond with jobHash after successful upload
    res.status(200).send({ message: 'Files uploaded successfully.', jobHash });
  } catch (err) {
    console.error('Error during file upload:', err);
    // Remove erroneous job from uploads
    if (jobHash) {
      if (fs.existsSync(jobUploadDir)) {
        fs.rmSync(jobUploadDir, { recursive: true });
      }
    }
    res.status(500).send('Internal Server Error');
  }
}
