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

const UPLOAD_DIR = path.resolve('./uploads');
const JOB_DIR = path.resolve('./public/jobs');

// Data structure to store pending uploads
const pendingUploads = {};

export default async function uploadHandler(req, res) {
  try {
    if (req.method !== 'POST') {
      return res.status(405).end('Method Not Allowed');
    }

    // Generate a unique hash for the job
    let jobHash;
    let jobUploadDir;
    do {
      jobHash = crypto.randomBytes(8).toString('hex');
      jobUploadDir = path.join(UPLOAD_DIR, jobHash);
    } while (fs.existsSync(jobUploadDir));

    // Add job to pending uploads
    pendingUploads[jobHash] = true;

    // Create job specific directories
    const genomeDir = path.join(jobUploadDir, 'genome');
    const gffDir = path.join(jobUploadDir, 'gff');

    try {
      await fs.promises.mkdir(jobUploadDir, { recursive: true });
      await fs.promises.mkdir(genomeDir, { recursive: true });
      await fs.promises.mkdir(gffDir, { recursive: true });
    } catch (err) {
      console.error('Error during directory creation:', err);
      throw err;
    }

    const form = new IncomingForm({ hashAlgorithm: 'sha256' });
    form.uploadDir = jobUploadDir;
    form.keepExtensions = true;
    form.multiples = true;
    form.maxFileSize = MAX_FILE_SIZE;
    form.maxTotalFileSize = MAX_TOTAL_FILE_SIZE;

    let masterTableFile;
    let motifNumber;
    let clientHashes;
    let uploadError;

    const uploadTimeout = calculateUploadTimeout(); // in milliseconds

    const formPromise = new Promise((resolve, reject) => {
      const timer = setTimeout(() => {
        uploadError = new Error('Upload timeout exceeded');
        form.emit('abort');
        reject(uploadError);
      }, uploadTimeout);

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

      form.parse(req, (err, fields, files) => {
        clearTimeout(timer); // Clear the timeout if the form parsing finishes in time
        if (err) {
          reject(err);
        } else if (uploadError) {
          reject(uploadError);
        } else {
          motifNumber = fields.motifNumber;
          clientHashes = JSON.parse(fields.fileHashes);
          // Compare the hashes
          for (let field in files) {
            for (let file of files[field]) {
              const formHash = file.hash;
              const clientHash = clientHashes[file.originalFilename];
                if (formHash !== clientHash) {
                  uploadError = new Error('Hash mismatch detected');
                  reject(uploadError);
                  return;
                }
            }
          }
          resolve({ fields, files });
        }
      });
    });

    try {
      await formPromise;
    } catch (err) {
      console.error('Error during form parsing:', err);
      throw err;
    }

    if (uploadError) {
      // Remove job from pending uploads
      delete pendingUploads[jobHash];
      throw uploadError;
    }

    // Create jobDir after successful upload
    const jobDir = path.join(JOB_DIR, jobHash);
    try {
      await fs.promises.mkdir(jobDir, { recursive: true });
    } catch (err) {
      console.error('Error during job directory creation:', err);
      throw err;
    }

    // Create params.json after successful upload
    const params = {
      inputTable: path.join(jobUploadDir, masterTableFile),
      inputGenomes: genomeDir,
      inputGFFs: gffDir,
      motifNumber: motifNumber,
      outputDir: jobDir,
    };
    try {
      await fs.promises.writeFile(path.join(jobDir, 'params.json'), JSON.stringify(params));
    } catch (err) {
      console.error('Error during params.json file creation:', err);
      throw err;
    }

    // Respond with jobHash after successful upload
    console.log(jobHash);
    res.status(200).send({ message: 'Files uploaded successfully.', jobHash: jobHash });

  } catch (err) {
    // Remove job from pending uploads
    delete pendingUploads[jobHash];
    console.error('Error during file upload:', err);
    res.status(500).send('Internal Server Error');
  }
}
