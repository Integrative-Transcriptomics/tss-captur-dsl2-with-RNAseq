import React, { useState } from 'react';
import Link from 'next/link';
import { useRouter } from 'next/router';
import styles from '../styles/Sidebar.module.css';

const Sidebar = ({ isOpen }) => {
  const [isLoggedIn, setIsLoggedIn] = useState(false);
  const router = useRouter();

  // Mockup login handlers
  const handleLogin = () => {
    setIsLoggedIn(true);
  };

  const handleLogout = () => {
    setIsLoggedIn(false);
    router.push('/');
  };

  return (
    <div className={`${styles.sidebar} ${isOpen ? styles.open : ''}`}>
      {isLoggedIn ? (
        // Logged in
        <nav className={styles.navLinks}>
          {/* Navigation links*/}
          <Link href="/dashboard">Dashboard</Link>
          <Link href="/newjob">New Job</Link>
          <button onClick={handleLogout}>Logout</button>
        </nav>
      ) : (
        // Not logged in
        <div className={styles.loginField}>
          <form>
            <div>
              <label htmlFor="email">Email:</label>
              <input type="email" id="email" name="email" required />
            </div>
            <div>
              <label htmlFor="password">Password:</label>
              <input type="password" id="password" name="password" required />
            </div>
            <div>
              <button type="submit" onClick={handleLogin}>Login</button>
            </div>
            <div>
              {/* Navigation links*/}
              <Link href="/register">Register</Link>
            </div>
          </form>
        </div>
      )}
    </div>
  );
};

export default Sidebar;