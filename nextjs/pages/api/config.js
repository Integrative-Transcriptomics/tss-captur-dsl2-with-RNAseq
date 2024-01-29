// Define upload configuration
const MAX_FILE_SIZE = (200) * 1024 * 1024; // 200MB
const MAX_TOTAL_FILE_SIZE = (500) * 1024 * 1024; // 500MB

export default function configHandler(req, res) {
  try {
    if (req.method === 'GET') {
      res.status(200).json({ maxFileSize: MAX_FILE_SIZE, maxTotalFileSize: MAX_TOTAL_FILE_SIZE });
    } else {
      res.status(405).json({ error: 'Method not allowed' });
    }
  } catch (error) {
    res.status(500).json({ error: 'An unexpected error occurred' });
  }
}

