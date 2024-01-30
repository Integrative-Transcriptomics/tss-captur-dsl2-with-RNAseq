import Link from 'next/link';

const SidebarLinks = () => {
  return (
    <div className="sidebar-links">
      <ul className="navbar-nav">
        {/* File upload */}
        <li className="nav-item">
          <Link className="nav-link" href="/">
              <i className="fas fa-home"></i>
              <span>File Upload</span>
          </Link>
        </li>
        {/* Example data */}
        <li className="nav-item">
          <Link className="nav-link" href="/example">
              <i className="fas fa-database"></i>
              <span>Example Data</span>
          </Link>
        </li>        
        {/* About page */}
        <li className="nav-item">
          <Link className="nav-link" href="/about">
              <i className="fas fa-info-circle"></i>
              <span>About</span>
          </Link>
        </li>
      </ul>
    </div>
  );
};

export default SidebarLinks;