import Layout from '../components/mockup/Layout';
import Link from 'next/link';

export default function DashboardPage() {
  // Mock job data
  const jobs = [
    { id: 1, name: 'Job 1', progress: 'In Progress', resultsLink: '/results/1' },
    { id: 2, name: 'Job 2', progress: 'Completed', resultsLink: '/results/2' },
    { id: 3, name: 'Job 3', progress: 'Not Started', resultsLink: '/results/3' },
  ];

  return (
    <Layout>
      <div>
        <h1>Dashboard</h1>
        {/* Display jobs */}
        {jobs.map((job) => (
          <div key={job.id}>
            <h2>{job.name}</h2>
            <p>{job.progress}</p>
            {job.progress === 'Completed' && <Link href={job.resultsLink}>View Results</Link>}
          </div>
        ))}
      </div>
    </Layout>
  );
}