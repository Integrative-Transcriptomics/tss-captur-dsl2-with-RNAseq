import React from 'react';
import Head from 'next/head';
import Link from 'next/link';
import Script from 'next/script';
import SidebarLinks from './SidebarLinks';

const Layout = ({ children }) => {
  return (
    <>
      <Head>
        <title>TSS-Captur</title>
        <meta charset="utf-8" />
        <meta http-equiv="X-UA-Compatible" content="IE=edge" />
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no" />
        <meta name="description" content="" />
        <meta name="author" content="" />
     </Head>
       <div id="wrapper">
        {/* Sidebar */}
        <ul className="navbar-nav bg-gradient-primary sidebar sidebar-dark accordion" id="accordionSidebar">
          {/* Sidebar - Brand */}
          <Link className="sidebar-brand d-flex align-items-center justify-content-center" href="/">
              <div className="sidebar-brand-icon rotate-n-15">
                <i className="fas fa-dna"></i>
              </div>
              <div className="sidebar-brand-text mx-3">TSS-Captur</div>
          </Link>
          <SidebarLinks />
          {/* Divider */}
          <hr className="sidebar-divider my-0" />
          <hr className="sidebar-divider" />
          <div className="text-center d-none d-md-inline">
              <button className="rounded-circle border-0" id="sidebarToggle"></button>
          </div>
        </ul>
        {/* Content Wrapper */}
        <div id="content-wrapper" className="d-flex flex-column">
          {/* Main Content */}
          <div id="content">
            <div className="container-fluid">
              {children}
            </div>
          </div>
          {/* Footer */}
          <footer className="sticky-footer bg-white">
            <div className="container my-auto">
              <div className="copyright text-center my-auto">
                <span>Copyright &copy; TSS-CAPTUR 2021</span>
              </div>
            </div>
          </footer>
        </div>
      </div>
      {/* Custom scripts for all pages*/}
      <Script src="/static/vendor/jquery/jquery.min.js" strategy="beforeInteractive"></Script>
      <Script src="/static/vendor/bootstrap/js/bootstrap.bundle.min.js" strategy="beforeInteractive"></Script>
      <Script src="/static/vendor/jquery-easing/jquery.easing.min.js" strategy="beforeInteractive"></Script>
      <Script src="/static/js/sb-admin-2.min.js" strategy="afterInteractive"></Script>
    </>
  );
};

export default Layout;