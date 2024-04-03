import React from 'react';
import Head from 'next/head';
import Link from 'next/link';
import Script from 'next/script';
import SidebarLinks from './SidebarLinks';

// Main layout component
const Layout = ({ children }) => {
  return (
    <>
      {/* Head element for metadata */}
      <Head>
        <title>TSS-Captur</title>
        <meta charSet="utf-8" />
        <meta httpEquiv="X-UA-Compatible" content="IE=edge" />
        <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no" />
        <meta name="description" content="" />
        <meta name="author" content="" />
      </Head>
      <div id="wrapper">
        {/* Sidebar navigation */}
        <ul className="navbar-nav bg-gradient-primary sidebar sidebar-dark accordion" id="accordionSidebar">
          {/* Brand logo and name */}
          <Link className="sidebar-brand d-flex align-items-center justify-content-center" href="/">
            <div className="sidebar-brand-icon rotate-n-15">
              <i className="fas fa-dna"></i>
            </div>
            <div className="sidebar-brand-text mx-3">TSS-Captur</div>
          </Link>
          {/* Sidebar links component */}
          <SidebarLinks />
          {/* Dividers for separating sections */}
          <hr className="sidebar-divider my-0" />
          <hr className="sidebar-divider" />
          {/* Sidebar toggle button for mobile view */}
          <div className="text-center d-none d-md-inline">
            <button className="rounded-circle border-0" id="sidebarToggle"></button>
          </div>
        </ul>
        {/* Main content area */}
        <div id="content-wrapper" className="d-flex flex-column">
          {/* Container for child components */}
          <div id="content" style={{ maxWidth: '75%', marginTop: '20px' }}>
            <div className="container-fluid d-flex flex-column">
              {children}
            </div>
          </div>
          {/* Footer*/}
          <footer className="sticky-footer bg-white">
            <div className="container my-auto">
              <div className="copyright text-center my-auto">
                <span>University of Tübingen – <a href="https://uni-tuebingen.de/fakultaeten/mathematisch-naturwissenschaftliche-fakultaet/fachbereiche/informatik/lehrstuehle/integrative-transkriptomik/home/">Integrative Transcriptomics <i className="fas fa-external-link-alt"></i></a> –TSS-CAPTUR 2024</span>


              </div>
            </div>
          </footer>
        </div>
      </div>
      {/* Scripts for functionality and interactivity */}
      <Script src="/static/vendor/jquery/jquery.min.js" strategy="beforeInteractive"></Script>
      <Script src="/static/vendor/bootstrap/js/bootstrap.bundle.min.js" strategy="beforeInteractive"></Script>
      <Script src="/static/vendor/jquery-easing/jquery.easing.min.js" strategy="beforeInteractive"></Script>
      <Script src="/static/js/sb-admin-2.min.js" strategy="afterInteractive"></Script>
    </>
  );
};

export default Layout;