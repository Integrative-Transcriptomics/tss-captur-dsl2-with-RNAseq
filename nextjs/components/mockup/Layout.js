import React, { useState } from 'react';
import Sidebar from './Sidebar';

const Layout = ({ children }) => {
  const [isOpen, setIsOpen] = useState(false);
  const toggleSidebar = () => setIsOpen(!isOpen);

  return (
    <>
      <button style={{position: 'absolute', top: 0, left: isOpen ? '200px' : 0}} onClick={toggleSidebar}>
        Toggle Sidebar
      </button>
      <Sidebar isOpen={isOpen} />
      <main>{children}</main>
    </>
  );
};

export default Layout;