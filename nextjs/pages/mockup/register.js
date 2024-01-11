import React, { useState } from 'react';
import Layout from '../components/mockup/Layout';

export default function AuthPage() {
  const [email, setEmail] = useState('');
  const [password, setPassword] = useState('');
  const [reEnterEmail, setReEnterEmail] = useState('');
  const [reEnterPassword, setReEnterPassword] = useState('');
  const [emailError, setEmailError] = useState('');
  const [passwordError, setPasswordError] = useState('');
  const [reEnterEmailError, setReEnterEmailError] = useState('');
  const [reEnterPasswordError, setReEnterPasswordError] = useState('');

  const validateEmail = (email) => {
    const re = /\S+@\S+\.\S+/;
    return re.test(email);
  };

  const validatePassword = (password) => {
    return password.length >= 8;
  };

  const handleRegister = (e) => {
    e.preventDefault();
    setEmailError('');
    setPasswordError('');
    setReEnterEmailError('');
    setReEnterPasswordError('');

    if (!validateEmail(email)) {
      setEmailError('Invalid email format');
    }

    if (!validatePassword(password)) {
      setPasswordError('Password must be at least 8 characters long');
    }

    if (email !== reEnterEmail) {
      setReEnterEmailError('Emails do not match');
    }

    if (password !== reEnterPassword) {
      setReEnterPasswordError('Passwords do not match');
    }

    if (validateEmail(email) && validatePassword(password) && email === reEnterEmail && password === reEnterPassword) {
      // authentication
    }
  };

  return (
    <Layout>
      <h1>Register</h1>
      <form onSubmit={handleRegister}>
        <input type="email" value={email} onChange={(e) => setEmail(e.target.value)} placeholder="Email" required />
        {emailError && <p>{emailError}</p>}
        <input type="email" value={reEnterEmail} onChange={(e) => setReEnterEmail(e.target.value)} placeholder="Re-enter Email" required />
        {reEnterEmailError && <p>{reEnterEmailError}</p>}
        <input type="password" value={password} onChange={(e) => setPassword(e.target.value)} placeholder="Password" required />
        {passwordError && <p>{passwordError}</p>}
        <input type="password" value={reEnterPassword} onChange={(e) => setReEnterPassword(e.target.value)} placeholder="Re-enter Password" required />
        {reEnterPasswordError && <p>{reEnterPasswordError}</p>}
        <button type="submit">Register</button>
      </form>
    </Layout>
  );
}