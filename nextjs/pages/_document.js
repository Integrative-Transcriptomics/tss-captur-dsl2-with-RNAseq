import Document, { Html, Head, Main, NextScript } from 'next/document';

// Custom Document component to modify applications <Head> tags
class MyDocument extends Document {
  // Render the structure of the document
  render() {
    return (
      <Html>
        <Head>
          {/* Google Fonts Link for Nunito */}
          <link href="https://fonts.googleapis.com/css?family=Nunito:200,200i,300,300i,400,400i,600,600i,700,700i,800,800i,900,900i&display=swap" rel="stylesheet" />
          {/* External CSS Stylesheets */}
          <link href="/static/vendor/fontawesome-free/css/all.min.css" rel="stylesheet" type="text/css" />
          <link href="/static/css/sb-admin-2.min.css" rel="stylesheet" />
        </Head>
        <body>
          {/* Main content of the document */}
          <Main />
          <NextScript />
        </body>
      </Html>
    );
  }
}

export default MyDocument;
